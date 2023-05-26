# ROOPERT (Regenerativley cOOled Propulsion: Engine Research Toolbox)
Analysis for UW Society for Advanced Rocket Propulsion engine development, geared towards development of Liquid Courage 3.

STEP BY STEP FOR GETTING THIS TO RUN
1. You need to run this on linux. If you already have a linux boot, skip this step. Otherwise, I would recommend installing WSL. (https://learn.microsoft.com/en-us/windows/wsl/install). This is on windows. On max, there are other ways to get a linux kernel running which you will have to figure out with some googling. If neither of these options work (I would reallllly try to get WSL working), install a linux virtual machine. go to 
https://www.howtogeek.com/170870/5-ways-to-run-linux-software-on-windows/#:~:text=Virtual%20machines%20allow%20you%20to,it%20on%20a%20standard%20computer.
and follow the instructions. I use the vmware virtual machine, but you can install a seperate boot if thats easier. I HIGHLY RECOMMEND USING VMWARE if youre using a virtual machine.
2. Install an IDE on the linux VM or boot, like VScode or pycharm (I prefer VScode). If you're on WSL, just install vscode on your machine in windows. Once you get the IDE up, just look up "integrate WSL in ** IDE". For VScode, you need to install the WSL package (https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl). IDEs are like fancy text editors that make it wayyy easier to write and run code, as well as integrate GIT.
3. Once you get your IDE up, there should be away to clone a repository from git (version control) on the starting screen.
Make sure to seelect the "clone from URL" or "clone from HTTPS" (not the SSH) option. You should also make sure you are in a WSL or linux terminal (https://code.visualstudio.com/docs/remote/wsl) first.
Cloning a repository means you make a copy of the main directory on your machine. As you make changes to that directory,
it will only be on your personal machine and won't affect other people's work. If you want to send your changes to the
group you have to commit and push your code to the branch you're on. Then you can merge that branch with the main branch
to the main branch so everyone is working off your changes. GIT is really cool so google and learn more.
4. Go to Gitlab on your browser and navigate to the ROOPERT project. Then click the repository link on the left side of the screen.
Go to the top right where it says CLONE in a blue blox and click that button. Copy the HTTPS link. Paste that into wherever
it asks for the link on your IDE. HTTPs is way easier! It will prompt you for your username and password.
5. In order for CEA functions to run you need to isntall rocketCEA and rocketProps packages, 
which I could not get working on windows. I would recommend using ubuntu and following the directions 
at https://rocketcea.readthedocs.io/en/latest/quickstart.html. You need a foertran compiler, so just do what it says. Make sure you do it in a WSL terminal (you can open one from vscode or however else). Once you get the fortran compiler
set up you can just pip install both packages!
6. Python uses packages to do cool things other than just add numbers. You need those packages installed in the virtual environment
that your IDE is running your directory in. To make sure you're installing in the virtual env your IDE is running,
go to the terminal in your IDE. After installing the fortran compiler above,
To quickly install all required packages, once in the ROOPERT directory and inthe terminal (not the python console),
just run pip install -r requirements.txt
6. Once you get the repository cloned, use the git interface on your IDE to make a new branch so that you're not just working
on the main branch and can do whatever you want without breaking others' work. Pull from the main branch often to make sure your
code is up to date with what others are working on, and push to main when you get your stuff working.
7. Currently the code is very disorganized, but OoeyGoey.py will pop up the GUI. Other than that, Main.py has a lot of the stuff in it.

8. If you can't run either of these because of the error "line 5462, end of file", this is an error in the transport property calcluation from rocketCEA.
IDK why it was broken but it is in the new version so just install the previous version of rocketCEA (1.1.28) using pip.
If you get the error "module Toolbox.constant can't be found" or some similar error where python can't find the folders, look at the top of main to see what you need to add to your path to get python to see the files.


To quickly install all required packages, once in the ROOPERT directory, just run pip install -r requirements.txt

Planned development:

Phase 1: First order calcs working with mass approx to give general paramters of a valid rocket for other teams to begin iterating and giving us more specific numbers

Phase 2:  
- Develop low-level thrust chamber analysis, including cooling passage geometry, thrust chamber geometry, flow inside chamber, injector parameters, coupling of cooling, injection, and performance, etc. 
- Work closely with solidworks to make sure parameters are quickly passable to create a geometry (probably need seperate stuff for brazed tubng, milled chanels, coaxial shell, and whatever else). 
- Concurrently begin to set up skeleton code and reasonable interface for toher teams to put in their own analysis (RF for tanks and filling, structures for mass approximations, system interconnects for valve sizing and pressure drop, whoever else wants to contribute to this mess)
- find a reasonable way to store iterations and data in a compact and useable format.

Phase 3: Iterate design and keep adding complexity and resolution until you can just press a button and you get a rocket


