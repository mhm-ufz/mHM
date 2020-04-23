# mHM check cases

With the provided python script `run_mhm_checks.py` you can run all checks in the `case_**` folders with one or more given mhm executable(s).

It will be checked, if mhm finishes normally and if the output matches the reference values.

A summary for all cases is given at the end.

## Usage

    python run_mhm_checks.py [-h] [-e EXE [EXE ...]] [-v] [-l LOG_PATH] [-m MPI_NOP] [-t OPENMP_THREADS]

Run the mhm check cases with a given mhm executable.

## Optional Arguments

    -h, --help            show this help message and exit
    -e EXE [EXE ...], --exe EXE [EXE ...]
                          Paths to mhm exe[s]. (default: ['../mhm'])
    -v, --verbose         Show the mhm output. (default: no)
    -l LOG_PATH, --log_path LOG_PATH
                          Directory for mhm-logs. (default: the resp. case dir)
    -m MPI_NOP, --mpi MPI_NOP
                          Number of mpi processes. No openMP allowed! (default:
                          0)
    -t OPENMP_THREADS, --threads OPENMP_THREADS
                          Number of threads for openMP. No mpi allowed!
                          (default: 0)
    -s [SKIP [SKIP ...]], --skip [SKIP [SKIP ...]]
                          skip cases (case_01 case_03 ..) (default: [])

## Examples
Run mhm from parent directory in verbosity mode with mpi on 4 processes:

        python run_mhm_checks.py -e ../mhm -v -m 4

Silently run mhm (given with an absolute path) with openmp on 4 threads:

        python run_mhm_checks.py -e /abspath/mhm_openmp -t 4

Silently run mhm from parent directory:

        python run_mhm_checks.py

Run with multiple mhm exes:

        python run_mhm_checks.py -e ../mhm1 ../mhm2

## Cleanup
To remove the created output of mHM run:

        python clean_mhm_checks.py

## Author
    Sebastian Mueller

## Contributors
    Stephan Thober, Robert Schweppe

Written Feb. 2020.