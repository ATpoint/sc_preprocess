# containerlog

Build container with `docker build -t user/reponame:tag .` given the Dockerfile is in the current directory. Then, given the repository was created at DockerHub, push it to the Hub via `docker push user/reponame:tag`. This can either be done from local machine or directly via GitPod. In the future we might add a GitHub Actions to automate that upon pushes to the `environment.yml`.

## v1.2.0
- switched to micromamba base image

## v1.1.1
- updated to Bioc-3.13 and R-4.1

## v1.0.0
- first version 