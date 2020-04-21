from os import system as os
import sys

#Clear data directory


#Specify compilation settings
settings = ' -ffree-form -ffree-line-length-0 -fdefault-real-8 -O3 -Wunused-label '

#Compile all modules
os ("rm -rf mod/*")
os ("rm /Users/tomkimpson/Data/ThesisData/RT/*.txt")
os("gfortran -J mod/ -c"+settings+"parameters.f -o mod/1.o") 
os("gfortran -J mod/ -c"+settings+"constants.f -o mod/2.o") 
os("gfortran -J mod/ -c"+settings+"ODEs.f -o mod/3.o") 
os("gfortran -J mod/ -c"+settings+"initial_conditions.f -o mod/4.o") 
os("gfortran -J mod/ -c"+settings+"numerical_methods.f -o mod/5.o") 
os("gfortran -J mod/ -c"+settings+"IO.f -o mod/6.o") 
os("gfortran -J mod/ -c"+settings+"ray.f -o mod/7.o") 
os("gfortran -J mod/ -c"+settings+"optimization.f -o mod/8.o") 
os("gfortran -J mod/ -c"+settings+"main.f -o mod/9.o") 



#Link all modules
print ('Starting compilation')
os("gfortran mod/*.o -o GO")
print ('Compiled.') 


#Run the code
os("./GO")


