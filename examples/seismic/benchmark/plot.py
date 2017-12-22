import os
import sys
from subprocess import call
import re

#
# Expected directory tree:
#
# sessionname (eg SISC-paper-2018; to be provided by user)
#  | -- problemname (eg: acoustic, tti)
#  |       | -- backend (eg: core, yask)
#  |       |      | -- benchmode (eg: bmO2, bmdse)
#  |       |      |       | -- archname (eg: knl7250)
#  |       |      |       |       | -- grid512 -- files
#  |       |      |       |       | -- grid768 -- files
#  |       |      |       |       | -- grid1024 -- files
#
# E.g.: paper-SISC-2018/acoustic/core/bmO2/knl7250/grid512
#

if len(sys.argv) != 2:
    print("Usage: python plot.py absolute/path/to/session/name")

machines = {}
exec(open("machines.txt").read(), {}, machines)

subfolders = ['problemname', 'backend', 'benchmode', 'archname', 'grid']

benchmarkpy = os.path.join(os.environ['DEVITO_HOME'], 'examples', 'seismic',
                           'benchmark', 'benchmark.py')

root = sys.argv[1]
dirs = [x[0] for x in os.walk(root)]

root_depth = len(root.split('/'))
plot_dirs = [i for i in dirs if len(i.split('/')) - len(subfolders) == root_depth]

for i in plot_dirs:
    problem, backend, benchmode, arch, _ = i.split('/')[root_depth:]
    files = [f for f in os.listdir(i) if f.endswith('json')]
    if not files:
        print("WTF?")
        sys.exit(0)

    if arch not in machines:
        print("Don't own any information for %s architecture" % arch)
        sys.exit(0)
    machine = machines[arch]

    grid = re.search(r"shape\[([A-Za-z0-9_,]+)\]", files[0]).group(1).split(',')
    space_orders = sorted(set([int(re.search(r"so\[(\w+)\]", j).group(1)) for j in files]))
    tn = re.search(r"tn\[(\w+)\]", files[0]).group(1)

    os.environ['DEVITO_ARCH'] = 'intel'
    os.environ['DEVITO_BACKEND'] = backend
    os.environ['DEVITO_AUTOTUNING'] = 'aggressive'

    args = ['python', benchmarkpy, 'plot', '-P', problem, '-bm', benchmode]  # plot cmd
    args += ['-d'] + grid  # grid shape
    for j in space_orders:
        args += ['-so', str(j)]  # space orders tested
    args += ['-to', '2', '--tn', tn]  # time-related simulation parameters
    args += ['-r', i]  # plot dir
    args += ['--arch', arch, '--max-bw', str(machine['dram-stream-bw']),  # arch params
             '--flop-ceil', str(machine['machine-peak']), 'ideal',
             '--flop-ceil', str(machine['linpack-peak']), 'linpack']
    call(' '.join(args), shell=True)
