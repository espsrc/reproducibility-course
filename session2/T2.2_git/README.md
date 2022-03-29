# Introduction to control versions with `git`

:bulb: Objectives:

- Knowing what git is and the basic concepts of use.
- Manage the lifecycle and work with repositories under git.
- Operate with files and versions in the repository from git.
- Track changes to code using git.
- Interact and collaborate with other repositories using git.

<HR>

<!-- vscode-markdown-toc -->
* 1. [Content](#Content)
	* 1.1. [Terminology and operations](#Terminologyandoperations)
	* 1.2. [Creating and using a existing repository](#Creatingandusingaexistingrepository)
	* 1.3. [Adding content and modifiying a git repository](#Addingcontentandmodifiyingagitrepository)
	* 1.4. [Tracking changes](#Trackingchanges)
	* 1.5. [Basic workflow with ´git´](#Basicworkflowwithgit)
	* 1.6. [Working with branches](#Workingwithbranches)
	* 1.7. [Collaborating with other repositories](#Collaboratingwithotherrepositories)
* 2. [Exercises](#Exercises)
	* 2.1. [Creating a new repository and adding content](#Creatinganewrepositoryandaddingcontent)
	* 2.2. [Forking a repository to collaborate with it](#Forkingarepositorytocollaboratewithit)
	* 2.3. [Working with branches](#Workingwithbranches-1)

<!-- vscode-markdown-toc-config
	numbering=true
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

<HR>


##  1. <a name='Content'></a>Content

### Introduction to `git`

In this section we'll look at the basics of Git. We'll cover all the terminology and commands necessary for you to become a regular Git user. 

First, we need to explain what version control is. Version management can be seen and used as a way to keep track of changes to files in a clearly defined way. More importantly, we can access the history of these changes and get details of the modifications made to files over time. 

Modern version control systems also offer tools that allow code sharing, which makes working collaboratively easy and efficient.

The following examples do not correspond to a version management system.

<image #1 and image #2>

- false sense of security: it's easy to delete important things
- easy to make changes in the wrong directory
- not suitable for collaborations
- impossible to find the right version of the algorithm you wrote 3 months ago
- there are tools that make it easy
- tethered to a physical system and no decentralised or cloud support

Git is one of many version control systems available.

- CVS, Mercurial, SVN, and other.

**So why `git`**

- works on Linux, Windows and macOS
- easy to install: all major Linux distributions have relatively new versions available (you can always compile from source - if you need the latest version), .exe installer for Windows and an installer for macOS
- uses the command-line interface (CLI), but graphical interfaces are available
- only a few commands are needed to get started
- it is very difficult to delete data permanently by mistake, you have to be very specific
- works both locally and remotely


###  1.1. <a name='Terminologyandoperations'></a>Terminology and operations

These are the 4 basic concepts you need to know to use `git`.

**Repository**

[git] repository/repo - what this can mean can be different to what it actually means. Repository is the .git directory inside your working folder, which contains all the required repository files. Quite often though, when someone refers to 'a repository' they mean the whole folder which contains the files you are working on.

**Stages of your repository**

- *modified* - you made changes to some of your files, you saved them with your editor of choice (e.g. vim), but Git is not yet aware of them. It detects that something changed, but has not yet stored it in its long-term memory (a database)
- *staged* - you tell Git which changes should be stored in the database next (you stage these files)
- *committed* - Git takes all the staged files you ordered it to keep track of and stores their snapshots (how they look at that point in time) into the database

**Branches**

branch - most of the time, code development is not linear: you change some things in one place, add a feature in another and maybe remove something in a completely different file. What if you have multiple collaborators working on separate features? In such a situation, they can all work on their own, personal versions of the code (branches), not modifying everyone else's contributions and possibly the production-ready code.

**Main or Master branch**

main/master branch - this is the branch you want the world to see and use. It should only include files that you are confident other people will be able to use and use them without major issues. It should contain production-ready and working code only - if you need to add and test a new feature, create a separate branch for it.

###  1.2. <a name='Creatingandusingaexistingrepository'></a>Creating and using a existing repository

There are two ways to work with git repositories: 

- :one: Start an empty repository and add files for version control
- :two: Work with an available repository

So you have a new codebase or a document that you have just started working on and decided to start using Git to keep track of all the changes you make. As we have already mentioned, it is easy to start using Git and you need only a single command to start a repository for your project (yes, it's really that easy). Go to your project's directory and run:

```
git init
Initialised empty Git repository ...
```

If you have a link to the GitHub repository, for example your coworker shared a link with you to some amazing plotting library and wants you to have a look at the source code, you can download the contents of this repository with a single command:

```
git clone <repository link>
Cloning into...
...
```

For example, to clone an  Git repository using `git` command and hosted on GitHub, you should run: 

*:bulb: Go to your repository website and copy the code of the URL of it by clicking the button `Code`*
```
$ git clone https://github.com/<your_user_name>/reproducibidad-csic.git
Cloning into 'reproducibidad-csic'...
...
```

This command downloads the repository and automatically places it inside a git folder (and it initializes your repository to use with `git`). 

:bulb: In the following steps we will interact with the repository through the `git` command.

###  1.3. <a name='Addingcontentandmodifiyingagitrepository'></a>Adding content and modifiying a git repository

To add content you will need to have a text editing tool (vi, vim, gedit, sublime text, textmate, bloc de notas, ...).

First we are going to create a file with the name `functions.py` where we will add the following text:

```
def datatext():
    print ("...")
```

Once this is done we save the file inside the folder of the previous repository that we have cloned and we go back to the terminal window, where we will use `git`. 

Now let's check what changes we have pending, to do this use `git`:

```
git status
```

Here you will see the new file we have added.

```
...
```

Now open the text editor again and select the `README.md` file of the repository we have cloned, and add the following text at the end of it:

Now we consult again from the terminal screen the status of the repository:

```
...
```

We can see that two changes appear:
1.- The new addition of the previous file created (`functions.py`)
2.- A modification of `README.md`.

After this we can fix the changes by adding the files or directories we want to be part of this commit. 

```
git add README.md
git add functions.py
```

And finally we create the commit of this group of files. 
This change is however not permanently added to the repository history yet! To do that, we have to explicitly *commit* the changes:

```
git commit -m "Updated function file and README"
```

:bulb: *The commit message is a short description of the changes that have been made to the repository since the last commit. In our case, we have added a single, very short file, so an equally short commit message "Added a hello.txt file" will suffice. However, if you are working on a larger project, possibly with multiple other developers and researchers, you will want to make your commit messages as informative as possible, to make sure everyone can understand the changes made to the files without the need to dig through the source code unnecessarily. Even if you are working on your own, you will forget all the details after some time. It is always a good practice to include a short title for a commit and then underneath it to write a short paragraph detailing all the changes made in this commit.*

Type the following to see the status again:

```
git status
```

*Here you can see that we have a new status with our commited files. They are ready to be pushed to our cloud repository in GitHub.*

Once this is done we can upload the changes to the GitHub repository, because if we don't do this the changes will only live on your computer and not in the cloud (the GitHub service). To do this:

```
git push
```

Once this is done, the changes will be uploaded to GitHub and we can check that this is the case by accessing our repository website.

###  1.4. <a name='Trackingchanges'></a>Tracking changes 

What if we wanted to review the changelog we have made to our repository so far? Git has powerful tools for checking the version and version history of our files, so it's very easy to keep track of everything we've done during the development of our work.

:bulb: *Remember that in the previous trainiing session, we saw how you can check your repository versions visually.*

To start playing with the `git` logs, we use:

```
git log
```

The output of this command is very rich :smile: :

![practice example 19](./media/practice2.1.png) 

In this screen you can see the `author`, `date` and `commit message`.

If you hold ENTER or use the cursor to scroll down in the output text, you can scroll through the entire Git history, though it can take some time. git log may seem like a simple command, all it does at the end of the day is display your commit messages. However, it is extremely powerful and comes with a lot of options. Therefore, it is important that you familiarise yourself with the most common options or the ones that you think will be most useful in your day-to-day work.

TBC.

###  1.5. <a name='Basicworkflowwithgit'></a>Basic workflow with ´git´

TBC.

![practice example 19](./media/practice2.2.jpg) 


###  1.6. <a name='Workingwithbranches'></a>Working with branches

###  1.7. <a name='Collaboratingwithotherrepositories'></a>Collaborating with other repositories

##  2. <a name='Exercises'></a>Exercises

###  2.1. <a name='Creatinganewrepositoryandaddingcontent'></a>Creating a new repository and adding content

###  2.2. <a name='Forkingarepositorytocollaboratewithit'></a>Forking a repository to collaborate with it

###  2.3. <a name='Workingwithbranches-1'></a>Working with branches








