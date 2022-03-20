# Introduction to control versions with `git`

:bulb: Objectives:

- Knowing what git is and the basic concepts of use.
- Manage the lifecycle and work with repositories under git.
- Operate with files and versions in the repository from git.
- Track changes to code using git.
- Interact and collaborate with other repositories using git.

## Content

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


### Terminology and operations

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

### Creating and using a existing repository

There are two ways to work with git repositories: 

- Start an empty repository and add files for version control
- Work with an available repository

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

For example, to clone an  Git repository using Git of course and hosted on GitHub, you should run

```
$ git clone https://github.com/manuparra/biometrics.git
Cloning into 'biometrics'...
...
```

This command downloads the repository and automatically places it inside a git folder (and it initializes your repository to use with `git`).

### Adding content and modifiying a git repository

### Tracking changes 

### Basic workflow with ´git´

### Working with branches

### Collaborating with other repositories

## Exercises




