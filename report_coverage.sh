#!/bin/bash
CODECOV_TOKEN="b6e7dfd4-94be-4484-bceb-b3d0090ea612"
``
REPO_TOKEN="b6e7dfd4-94be-4484-bceb-b3d0090ea612" julia -e 'using Pkg; Pkg.add("Coverage");cd(Pkg.dir("TopologyPreprocessing")); using Coverage;  Codecov.submit_token(Codecov.process_folder())'
# REPO_TOKEN="b6e7dfd4-94be-4484-bceb-b3d0090ea612" julia -e 'using Pkg; import TopologyPreprocessing; joinpath(dirname(pathof(TopologyPreprocessing)), "..", paths...); Pgk.activate("."); using Coverage;  Codecov.submit_token(Codecov.process_folder())'
