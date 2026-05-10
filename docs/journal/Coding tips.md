# Vscode tips
Alt + j followed by Alt + o opens the julia REPL
Ctrl + Space - forces autocomplete (In vscode)
Shift + Enter - Run a cell denoted by ##
Shift + Delete - Delete all plots

# Julia REPL
* "]" opens pkg "package manager"
* ";" opens bash
* "?" opens the helper

# Julia tips
* ! marks functions that alters the arguments/parameters/inputs
* && and || are both short circuited - meaning that if the first argument is false in the case of && the second argument is not evaluated, which may be useful in cases where the second argument is dependent on the first being true.
* Bools are numbers
* "*" concatenates strings
* $(Variable) in a string will concatinate it with the rest of the string
* arguments after ";" in function arguments are keywords
* "..." splats/unpacks a tuple
* @view gives a subArray which is a lazy copy, hence changes made are reflected in the array subarrayed
* display(plot)
* @edit function(var) to read the source code for a given type! Very nice for debugging
* @. will broadcast all function calls, add $ infront of function call to not broadcast

* Build documentation using julia --project=. docs/make.jl

# Array tips
* stack(vector) to create one big array from vector of matricies, eachslice(array) reverses this

# Linux tips
* Use screen -S <sessionname> to create a screen, which can run a process while disconnected from
server. Use ctrl+A followed by D to detach from the session before logging off. Use sreen -r <sessionname> 
to reconnect to the session and continue where left off. Usefull for overnight simulations

# Github tips
* To track a file ignored by .gitignore do git add -f \<file>
* GitHub actions with "on: workflow_dispatch:" allows the user to run it manually in the repos
Actions tab
* The CompatHelper.yml needs the user to go into Settings > Actions > General > **Allow GitHub Actions to create and approve pull requests** (toggle it on)

# GPU tips
* use mapreduce to do both map and reduce operation in one go.

# Juliaup
* If you are struggling with a folder having a constant julia version check 'juliaup override status'

# Advectra package tips
* Make sure the broadcasting is done correctly, it is easy to mess up re-assign a variable 
instead of writing to it elementwise in-place.
* Having first order derivatives in either the linear or non-linear term does not seem to matter
* When Fourier transforming using dealiasing the dealiasing_coefficient needs to be divided by, while when inverse Fourier transforming the field needs to be multiplied by the dealiasing_coefficent.