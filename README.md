# GeoResistivity Processing Manual

## Installation

The programm is a distributed as a Python script. You need to have Python (3.4 or 2.7) installed on your computer (installed if you have ArcGIS).
There are dependencies to install separatelly:
- GDAL
- PyQt 4
- Numpy
- Matplotlib

They are common packages, check on Google how to have them installed on your machine.

## Usage
To run the script open the cmd terminal in windows and type <br/>
`"python /the/path/to/the/file.py"`<br/>
The windows will open, but you will stil have the terminal to find out about errors (nothing is perfect).


### Create new survey
To start with a new survey, you need to create a repertory with the following structure:

-SITENAME/ <br/>
    +output/ <br/>
        -geometry.txt <br/>
    +raw/<br/>

### Select the current site
In the space top left, lelect the directory of the site you work on wich contains the above mentionned folders. If it doesn't you will get an error message.

### Download data
To download data go to File>Import from RM85
  - select the number of the first grid you will download (it should be automatically the right one).
  - select the size of your grid (20 by 20)
  - select the port of the device
  - keep the baudrate at 9600
  - select the corner of the grid you used the most that day (you can still turn and flip everything around later, but it will save you time to have most of them all right)
  
###  Asemble grids

## Desctiption of lagorythms

### Despike

### Edge Matching

### Data previsualization

### Export
