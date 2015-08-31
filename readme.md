#Atomic Transposition of Atoms on 2D Surface Complexes 

This algorithm was implemented to transpose atoms oriented relative to a 2D atomic surface complex, set in a certain coordinate system, onto *another* 2D surface in a different coordinate system.

The algorithm uses Sympy and Math libraries.

#Data File Setup 

This script requires a data file set up in the following way, without spaces in between data.  File reading is performed by the script until a space is found, so any data that is relevant but not to be read can be stored after spaces, below the relevant data.  

- Surface 1, 3D Coordinate 1
- Surface 1, 3D Coordinate 2
- Surface 1, 3D Coordinate 3
- Surface 1, 3D Coordinate 4
- Surface 1, 3D Coordinate 5
- Surface 2, 3D Coordinate 1
- Surface 2, 3D Coordinate 2
- Surface 2, 3D Coordinate 3
- Surface 2, 3D Coordinate 4
- Surface 2, 3D Coordinate 5
- All atomic coordinates in terms of first coordinate system
	- Must be 4 columns, leading with atomic symbol 

** The first two coordinates from each surface must be RELATIVE TO EACH OTHER in order for transposition to be performed properly **

The script outputs 4 results, each corresponding to a different coordinate system relative to the original vector 

# How to Run the Script 
Script run with on the command line in the following way:

```python 
python transpose.py sampleData.txt 
```



