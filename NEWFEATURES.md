# New features 1909

## Several changes

- Folders take their names of the user input.
- Default objects are now numerated.
- Changes in default tube configuration.
- Adding dynamic library to work with Fedora.
- Some changes in the use of abs/rel paths to cross-platforms.
- Use labels of objects in Post Process.
- Change an object property on Modeling, reflect the change in Post Process.
- Generate default PP after run, but with question.
- Cleaner way to handle exceptions and show error messages.
- Check the size of the folders during the run, to refresh data in Post Process Tab in a more convinient way.
- Manage with assert empty data on plots.
- Change the ubication of the items in the basic wizard layout to fit the screen.
- Free plots only for cylinders and tanks.

## Corrected bugs

- Divide the initial state of the tank and cylinders of node zero and the rest.
- Adding cd ports with a tube - tank connection.
- Avoid malfunction if a valve key is missing.
- Different opening/closing angles on default valves.
- Don't depend on header.py. Activate the postPro at the start, only with the folder(s) located.
- Bug in x-conversion in free plots.
- Allow floats in tube length equispaced.
- Correction when erase an type of object with more than 10 elements.
- Invert the PV default diagram.
- Correction in units conversion from K to C.

## New features

- Added default cylinders in the wizard for 2/4 stroke and CI/SI.
- Inform simulated/to simulate rpms.
- Possibility to run and post-process at the same time.
- Completed ToolBar.
- Default Post Process Plots.
- Help Tab.
- Copy/Paste individual items.
- Save/Load particular configurations for items.
- Equispaced values on saving nodes for tubes.
- Function to change color/format/width of a particular curve.
- Edit curve attributes and delete curves using keyboard shortcuts.
- Progress bar for default Post Process.
