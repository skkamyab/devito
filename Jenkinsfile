pipeline {
  agent {
    docker { image 'ubuntu' }
  }
  environment {
    PATH = "/usr/local/bin:/usr/bin:/bin"
  }
  stages {
    stage('1') {
      steps {
        sh 'lsb_release -a'
      }
    }
  }
} 
