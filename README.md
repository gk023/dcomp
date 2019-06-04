###### Project directory ######

In this directory, the files needed to build and run the project are included.

More details are given in the documentation file about the specific functio of each file.

To build Boost:

   tar -xzvf boost_1_67_0.tar.gz

   Change makefiles LIBS variable to have -I ./boost_1_67_0. If you built boost somewhere
    outside this directory, just make sure the path after -I points to the top of the directory.

To build the code:

   make -f make_dcomp_main

To run the code: 

   python run_dcomp.py <path/to/card/file>. An example card file is given in ExampleCard.in

To plot the ouput:

   python plot.py <path/to/cross/section/files>. By default, the cross sections are saved in 
    cross_sections.out


g++ is the compiler I used to build the code. GSL will also need to be loaded to get the 
 integration routines
