# The XXI century lab book

In this section, we will see the value of Jupyter notebooks as a dynamic tool for exploratory analysis. We will learn how to initialize and navigate through notebooks, the basic structure and syntax to use a Jupyter notebook and the notebook cells.

## Introduction to Jupyter notebooks
The Jupyter Notebook is an open-source web application that allows you to create and share documents that contain code, equations, visualizations and text. The functionality is partly overlapping with [R Markdown](https://rmarkdown.rstudio.com/), in that they both use markdown and code chunks to generate reports that integrate results of computations with the code that generated them. Jupyter Notebook comes from the Python community while R Markdown was developed by RStudio, but you could use most common programming languages in either alternative. In practice though, it's quite common that R developers use Jupyter but probably not very common that Python developers use RStudio.

#### What are Jupyter notebooks for?
An excellent question! Some applications could be:

* Python is lacking a really good IDE for doing exploratory scientific data analysis, like RStudio or Matlab. Some people use it simply as an alternative for that.
* An early ambition with Jupyter notebooks was to be analogous to the lab notebook used in a wet lab. It would allow the scientist to document her day-to-day work and interweave results, ideas, and hypotheses with the code. From a reproducibility perspective, this is one of the main advantages.
* Jupyter notebooks can be used to provide a tighter connection between your data and your results by integrating the results of computations with the code that generated them. They can also do this in an interactive way that makes them very appealing for sharing with others.
* Notebooks are great tools for teaching. It is common now to share lecture materials for training schools as notebooks (see examples in the [Resources](#other-resources) section)
* The community around Jupyter notebooks is large and dynamic, and there are tons of tools for sharing, displaying or interacting with notebooks.
* In a research team, notebooks can be used to train new people on the methods and types of analysis used in your field.

As always, the best way is to try it out yourself and decide what to use it for!

#### Understanding the Jupyter nomenclature

 - [Project Jupyter](https://jupyter.org/): is the **project** that develops open-source software, open-standards, and services for interactive computing across dozens of programming languages.
 - A Jupyter notebook: The actual `.ipynb` **file format** that constitutes the notebooks. It is actually a json file.
 - Jupyter Notebook: is the **web application** that you use for creating, managing and running notebooks. Can be executed in a terminal with `jupyter notebook` [Try me live](https://mybinder.org/v2/gh/ipython/ipython-in-depth/master?filepath=binder/Index.ipynb)
 - [Jupyter Lab](https://jupyterlab.readthedocs.io/en/stable/index.html): is the next-generation web-based user interface for Project Jupyter. An advanced web application that can be used to work with notebooks but also text files, csv files, images, pdfs, terminals, etc. Can be executed in a terminal with `jupyter lab`. [Try me live](https://mybinder.org/v2/gh/jupyterlab/jupyterlab-demo/try.jupyter.org?urlpath=lab)
- [JupyterHub](https://jupyterhub.readthedocs.io/en/stable/): a multi-user version that can serve multiple instances of Jupyter Notebook servers to be used for example in a class of students, a corporate workgroup or a research lab.


## The Jupyter Notebook dashboard
One thing that sets Jupyter Notebook apart from what you might be used to is that it's a web application, i.e. you edit and run your code from your browser. But first you have to start the Jupyter Notebook server.

```no-highlight
$ jupyter notebook
[I 19:37:44.483 NotebookApp] Serving notebooks from local directory: /home/jmoldon/droplets_dev/droplets
[I 19:37:44.483 NotebookApp] The Jupyter Notebook is running at:
[I 19:37:44.483 NotebookApp] http://localhost:8888/?token=f024f1e7495d4c99125c860951f3b0a4679890e312937e54
[I 19:37:44.483 NotebookApp]  or http://127.0.0.1:8888/?token=f024f1e7495d4c99125c860951f3b0a4679890e312937e54
[I 19:37:44.483 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 19:37:44.519 NotebookApp]

    To access the notebook, open this file in a browser:
        file:///home/jmoldon/.local/share/jupyter/runtime/nbserver-5881-open.html
    Or copy and paste one of these URLs:
        http://localhost:8888/?token=f024f1e7495d4c99125c860951f3b0a4679890e312937e54
     or http://127.0.0.1:8888/?token=f024f1e7495d4c99125c860951f3b0a4679890e312937e54
```

Jupyter Notebook probably opened up a web browser for you automatically, otherwise go to the address specified in the message in the terminal. Note that the server is running locally (as [http://localhost:8888](http://localhost:8888)) so this does not require that you have an active internet connection. The notebook server will be initialized in the local directory it was executed on, so you will see the files in that directory. The working path is also displayed in the first line of the terminal. Everything you do in your Notebook session will be stored in this directory, so you won't lose any work if you shut down the server.

![](images/jupyter_dashboard.png)

What you're looking at is the Notebook dashboard. This is where you manage your files, notebooks, and kernels. The Files tab shows the files in your directory. The Running tab keeps track of all your processes. The third tab, Clusters, is used for parallel computing and won't be discussed further in this tutorial. You should see now the files of the droplets repository.

You can do several things: you can navigate to any file: for example, click on the file `environment.yml` to see its contents. You can close the browser tap when finished. You can also upload files or create new ones. You can even press New > Terminal to start a bash terminal in your web browser!

![](images/jupyter_file_terminal.gif)
## The very basics
Jupyter notebooks are made up out of cells, and you are currently standing in the first cell in your notebook. The fact that it has a green border indicates that it's in "Edit mode", so you can write stuff in it.
A blue border indicates "Command mode" (see below).
Cells in Jupyter notebooks can be of two types: markdown or code.

* **Markdown** - These cells contain static material such as captions, text, lists, images and so on. You express this using Markdown, which is a lightweight markup language. Markdown documents are plain text files that can then be converted to other formats for viewing (the document you're reading now is written in Markdown and then converted to HTML). More details in the next section [Some Markdown basics](#some-markdown-basics)

* **Code** - These are the cells that actually do something, just as code chunks do in R Markdown. You can write code in dozens of languages and all do all kinds of clever tricks. You then run the code cell and any output the code generates, such as text or figures, will be displayed beneath the cell. We will get back to this in much more detail, but for now it's enough to understand that code cells are for executing code that is interpreted by a kernel (in this case the Python version in your Conda environment).

Before we continue, here are some shortcuts that can be useful. Note that they are only applicable when in command mode (blue frames). Most of them are also available from the menus.
These shortcuts are also available from the **Help** menu in your notebook (there's even an option there to edit shortcuts).


* <kbd>Enter</kbd>: enter Edit mode
* <kbd>Esc</kbd>: enter Command mode
* <kbd>Ctrl</kbd>+<kbd>Enter</kbd>: run the cell
* <kbd>Shift</kbd>+<kbd>Enter</kbd>: run the cell and select the cell below
* <kbd>Alt</kbd>+<kbd>Enter</kbd>: run the cell and insert a new cell below
* <kbd>Ctrl</kbd>+s: save the notebook
* <kbd>Tab</kbd>: for code completion or indentation
* m/y: toggle between Markdown and Code cells
* dd: delete a cell
* a/b: insert cells above/below the current cell
* x/c/v: cut/copy/paste cells
* o: toggle output of the current cell

## Some Markdown basics

Markdown is easy to use, you can find details on the syntax here [Markdown syntax](https://www.markdownguide.org/basic-syntax/). A nice online resource to practise is [Dillinger](https://dillinger.io/)

Let's use our first cell to create a header. Change the format from
Code to Markdown in the drop-down list above the cell. Double click on
the cell to enter editing mode (green frame) and input "# My notebook"
("#" is used in Markdown for header 1). Run the cell with Shift-Enter.
Tada!

Markdown is a simple way to structure your notebook into sections with
descriptive notes, lists, links, images etc.

Below are some examples of what you can do in markdown. Paste all or parts
of it into one or more cells in your notebook to see how it renders. Make
sure you set the cell type to Markdown.

```
## Introduction
In this notebook I will try out some of the **fantastic** concepts of Jupyter Notebooks. I can even insert a link, follow [here](https://jupyter-notebook.readthedocs.io/en/stable/examples/Notebook/Working%20With%20Markdown%20Cells.html) 

### Markdown basics
Examples of text attributes are:

* *italics*
* **bold**
* `monospace`

I will need some latex here:

\begin{align}
\nabla \times \vec{\mathbf{B}} -\, \frac1c\, \frac{\partial\vec{\mathbf{E}}}{\partial t} & = \frac{4\pi}{c}\vec{\mathbf{j}} \\   \nabla \cdot \vec{\mathbf{E}} & = 4 \pi \rho \\
\nabla \times \vec{\mathbf{E}}\, +\, \frac1c\, \frac{\partial\vec{\mathbf{B}}}{\partial t} & = \vec{\mathbf{0}} \\
\nabla \cdot \vec{\mathbf{B}} & = 0
\end{align}

I will insert a very nice image here ![Moon](https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Lunar_libration_with_phase_Oct_2007_450px.gif/120px-Lunar_libration_with_phase_Oct_2007_450px.gif)

```

![](images/jupyter_markdown_cell.gif)


## Writing code
Now let's write some code! Since we chose a Python kernel, Python would be the native language to run in a cell. Enter this code in the second cell and run it:

```python
x = 2
print(x**2 + 3)
```

Note how the output is displayed below the cell. This interactive way of working is one of the things that sets Jupyter Notebook apart from normal python scripts. Python scripts are executed top-to-bottom in one run, while you work *in* a Jupyter notebook in a different way. You decide when to run a cell, so you are responsible for the order cells are executed. 

!!! warning
    Executing the cells in non-linear order is very common but also dangerous because the state of the program (the value of the variables and the functions) depend on the order you executed the cells. It is highly recommended that, from time to time, you restart the Kernel using "Restart & Clear Output" or "Restart & Run All" from the Kerner menu to make sure that the execution is linear.

Variables defined in cells become variables in the global namespace. Therefore, you can share information between cells. Try to define a function or variable in one cell and use it in the next. For example:

```python
import numpy as np
phase = np.linspace(0, 2*np.pi, 100)
```

and

```python
import matplotlib.pyplot as plt
plt.plot(phase, np.sin(phase))
```

Your notebook should now look something like this.

![](images/jupyter_basic_update.png)

The focus here is not on how to write Markdown or Python; you can make really pretty notebooks with Markdown and you can code whatever you want with Python. Rather, we will focus on the Jupyter Notebook features that allow you to do a little more than that.

## Additional notebook features

Here we show some features included in the notebooks to make your life easier. You can find this and more features [here](https://www.dataquest.io/blog/jupyter-notebook-tips-tricks-shortcuts/)

#### LaTeX formulas

As we have already seen, you can write formulas using the markdown cells. When you write LaTeX in a Markdown cell, it will be rendered as a formula using MathJax. Note that MathJax is a javascript program to render formulas but it will not have by default all the power of a full LaTeX installation.

This:

`$P(A \mid B) = \frac{P(B \mid A)P(A)}{P(B)}$`

becomes this:

$P(A \mid B) = \frac{P(B \mid A)P(A)}{P(B)}$

# Other resources

- Jake VanderPlas youtube series on Reproducible data analysis with jupyter [Youtube](https://www.youtube.com/playlist?list=PLYCpMb24GpOC704uO9svUrihl-HY1tTJJ)
- Try Jupyter in your browser [link](https://jupyter.org/try)
- Quickview Notebook sharing the Gravitational Wave detection [Notebook](https://github.com/losc-tutorial/quickview/blob/master/index.ipynb)
- A full Machine Learning course using Notebooks, for example [Lecture 1: Density Estimation](https://github.com/carmensg/IAA_School2019/blob/master/lectures/Day3-ZeljkoIvezic/notebooks/density_estimation.ipynb), [Lecture 3: Classification](https://github.com/carmensg/IAA_School2019/blob/master/lectures/Day3-ZeljkoIvezic/notebooks/classification.ipynb) and [Lecture 4: Dimensionality Reduction](https://github.com/carmensg/IAA_School2019/blob/master/lectures/Day3-ZeljkoIvezic/notebooks/dimensionality_reduction.ipynb).
- Another example of the full tutorial contents on an international Python conference: [PyCon 2015 Scikit-learn Tutorial](https://github.com/jakevdp/sklearn_pycon2015)

Not everything is a fairy tale in Jupyter notebooks. So when you are ready to not like Jupyter notebooks you can read these discussions:
- Blog post: [Why I Don't Like Jupyter](http://opiateforthemass.es/articles/why-i-dont-like-jupyter-fka-ipython-notebook/)
- Blog post: [Why Jupyter Is Not My Ideal Notebook](https://www.sicara.ai/blog/2019-02-25-why-jupyter-not-my-ideal-notebook)
- Youtube: [I don't like notebooks.- Joel Grus (Allen Institute for Artificial Intelligence)](https://www.youtube.com/watch?v=7jiPeIFXb6U)
