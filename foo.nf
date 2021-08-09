#! /usr(bin/env nextflow

nextflow.enable.dsl=2

def latest_sha = "git rev-parse HEAD".execute()

println latest_sha
