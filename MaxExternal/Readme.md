# Max MSP Explination:
The ladder filter model was built into various formats.

## Mac
The Max MSP external was built on top of the Min-Devkit framework. Refer to their GitHub for build instructions. Their repo must be cloned into the Packages directory. There are some extra configuration steps required to build the mac version:

1) Copy Osprey.ladder_tilde folder into the /source/projects directory from the location of the min-devkit clone (Typically /Documents/Max [v#]/Packages/min-devkit, although this could change if you created your own package as per the min-devkit instructions).
2) The CMakeLists.txt file might have to be slightly altered as per your Mac's configuraton. If you need a starting point, copy the CMake file from the min.xfade_tilde project.
3) signal_routing_objects.cpp and signal_routing_objects.h must still be present in the /source/projects/shared folder.
4) /source/min-api and /source/min-lib must still be present
5) Generate the XCode project as per the instructions on the Min-Devkit github (Mac instructions)
6) In XCode, change the macOS Deployment Target to 14.2 for both the project and for ALL_BUILD
7) Build with CMake as per the instructions on the Min-Devkit github
8) Files will be in the /externals folder

## Windows
The Windows build is simpler. Visual Studio Community and CMake are requirements.
1) Same as step 1 on Mac
2) Replace the CMakeLists.txt in the project folder with the one from the min.xfade_tilde project
3) Same as Mac
4) Same as Mac
5) Using CMake GUI, select source as the package folder (min-devkit, or whatever package you created). Select the /builds directory in that folder as the build path (create if necissary). Select the version of Visual Studio you are using, and generate the project.
6) Open the project in Visual Studio. Configure the builds to the relevant files and to be Release.
7) Build project in Visual Studio. Will be in /externals folder

## GenExpr
This doesn't need compiling and can be just be copied into your Max project directory. It can be opened in a Max project by creating an object called gen~ TransistorLadder. From here, it can by used or edited. Note that the CPU usage is much higher. There are also a couple minor differences:

1) The bias parameter messages are just called b
2) The feedback is not a signal input, just a parameter like business called fbk


