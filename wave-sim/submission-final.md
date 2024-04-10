## FEM (final submission)

Please fill this out and submit your work to Gradescope by the deadline.

### Output Video
![](final2.gif)

### Design Choices
- Extracting the surface mesh.
    
    I extracted the surface mesh by first creating a hashing function for faces that orders the face indices from least to greatest and turns them into a unique integer where the first number is multiplied by the total number of vertices squared, the second number is multiplied by the total number of vertices, and the largest number by 1 and then they're all added together. Then I created a hashmap from integers to faces. When I went to add a face in the mesh, I hashed the face, and checked if the hashmap already had the face; if it didn't I added the { hash, face } pair to the hashmap or if it did I deleted the face from the hashmap. Then when I had checked all the faces I iterated through all of the values of the hashmap and added them to an std::vector which I used to initialize the mesh. (Sorry I know that was more than 3 sentences, but I couldn't condense it).

- computing and applying internal forces
    I had a `Node` class that kept track of velocity, force, position, and acceleration; I had a `Face` class that kept track of its original area, normal, the nodes that define it and the node opposite it on the tetrahedron; finally, I had a `System` class that kept track of all the nodes and the tetrahedrons in the system. The tetrahedrons kept track of all of their faces, and the faces, again, pointed to the corresponding nodes, so that I could perform all the necessary calculations to calculate stress and strain in `updateCalculations` and actually apply them to the nodes in `updatePositions`.

- collision resolution
    I had a method in the `System` class called `resolveCollisions` which used implicit equations to check if the body was inside the collider.

- your explicit integration method
    I used the midpoint method by first saving the current positions and velocities, updating calculations, stepping forward by 1 half of the time step, updating calculations, resetting the positions and velocities, and stepping forward with `updatePositions` by one time step.
    I also did adaptive time step detailed below.

### Extra Features 
Briefly explain your implementation of any extra features, provide output images, and describe what each image demonstrates.

#### Make the visualizer pretty (5 points)

I implemented a shader that would show how forces were applied throughout the body at different tetrahedrons.

![](finalForce.gif)

#### Adaptive time stepping (5 points)

I implemented adaptive time stepping by saving the initial positions and velocities, stepping a whole step, storing all positions and velocities in a single vector, reverting, stepping forward two half steps, getting the state from all positions and velocities, and performing calculations accordingly according to the paper.

### Collaboration/References

I collaborated with classmates in collaborative hours including:

<ul>

<li> Autumn Tilley (atilley) </li>
<li> Jean Yoo </li>
<li> Jamie Chen </li>
<li> Luke Riley </li>
<li> Dylan Lee </li>

</ul>

### Known Bugs

None