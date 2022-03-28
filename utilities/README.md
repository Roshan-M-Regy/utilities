# CG simulation analysis codes
## Intrachain distances and polymer scaling exponent
```
python Chain_separation_s.py traj top(.data/.pdb) number_of_frames_to_skip
```

## Making movies in python
Requires the gif package
```
python plot_density_movies_mean_front.py 
```
This current example plots a 2x2 figure with contents changing in each window with a time ticker
Example ![Density-profile-gif](https://github.com/rmregy/analysis_codes/blob/main/example_movie.gif)

## Calculating oligomer sizes from gsd trajectory 
### Usage
```
python get_poly_index.py *.gsd 
```
### Inputs
Requires the following parameters to be set in script,<br/>
scanparam, parameter which was varied during scanning simulation<br/>
chainlen, length of chain (current script only works for a single component system)<br/>
headpos, positions of the head interface on a single chain<br/> 
tailpos, positions of the tail interface on a single chain<br/>
start, Which frame to start calculation from?<br/>
skip, How many frames to skip?<br/>

### Outputs
Writes a file "cluster_size_number.txt" containing 4 columns<br/>
Scanparam Largest_cluster Deviation number of clusters Deviation



# utilities
