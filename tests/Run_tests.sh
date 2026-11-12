#!/bin/bash
#
# Run the dabax tests on the code of this repository (not on an installed dabax).
# The download tests need network access (they are skipped otherwise).
#
# Usage:
#   ./tests/Run_tests.sh                   # uses "python" from the PATH
#   PYTHON=/path/to/python ./tests/Run_tests.sh
#
PYTHON=${PYTHON:-python}

# run from the repository root, so that "import dabax" loads the repository code
cd "$(dirname "$0")/.." || exit 1

echo "dabax tested: $($PYTHON -c 'import dabax.dabax_base as d; print(d.__file__)')"

$PYTHON -m unittest discover -s tests -v
# alternatively, with pytest:
# $PYTHON -m pytest tests -v
