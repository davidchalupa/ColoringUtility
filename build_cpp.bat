set PATH=C:\msys64\mingw64\bin;%PATH%
cmake -B build_cpp -G "MinGW Makefiles"
cmake --build build_cpp
