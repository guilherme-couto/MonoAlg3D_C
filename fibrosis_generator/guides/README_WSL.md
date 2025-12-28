WSL BUILD INSTRUCTIONS AND TROUBLESHOOTING
Target System: Windows Subsystem for Linux (WSL)
Dependencies: Octave, G++ (managed via Conda)

1. COMPILATION OF C++ SOURCES
The project contains C++ files (Octave2D.cpp, Octave3D.cpp) that must be compiled into MEX files callable by Octave. These files utilize the "mex.h" header (MATLAB API style).

Build Tool:
Use 'mkoctfile', which is provided by the Octave package in Conda.

Instructions:
1. Activate the Conda environment.
2. Navigate to the directory containing the .cpp files.
3. Run the following commands to generate the binary MEX files:

   $ mkoctfile --mex Octave2D.cpp
   $ mkoctfile --mex Octave3D.cpp

Expected Output:
This will generate files with extensions '.mex' (or '.mexa64'). These binaries are automatically detected by Octave as functions.

2. KNOWN ISSUE: OCTAVE TERMINAL FREEZE ON WSL
Description:
When running Octave via Conda on WSL, the standard 'octave' command may result in a frozen terminal. The prompt accepts no input, and text is not echoed back to the screen.

Cause:
This is a conflict between the Conda-provided 'readline' library and the system/WSL interface, often triggered by the GUI (Qt) initialization or terminal capability detection.

Workaround/Solution:
Do not run the standard 'octave' command. Instead, use the command-line interface (CLI) executable with specific flags to disable line editing features.

Command to run Octave interactively:
$ octave-cli --no-line-editing

Command to run scripts directly (Recommended):
$ octave-cli script_name.m

3. VERIFICATION
To verify the build was successful:
1. Open the terminal: $ octave-cli --no-line-editing
2. Check for MEX files: ls *.mex
3. Attempt to call the function helper (even if undocumented): help Octave2D

If the error returns "'Octave2D' is not documented" rather than "undefined," the binary has been successfully loaded.
