
# LODESTAR v1.0

## Building Executable

LODESTAR must be compiled with gcc, not clang.

### Linux

gcc is available on all linux distributions. LODESTAR is built by running make:

```
make
```

### MacOS

gcc can be installed on MacOS using homebrew. Once installed, run the following command to build LODESTAR:

```
make CC=$(ls /opt/homebrew/bin/ | grep '^gcc-..$')
```

## Rest of README Coming Soon!