# Shellfish Toxicity Modeling
Red tide is a type of harmful algae bloom that affects human and wild marine activities.
This project seeks to correlate the Hazardous Algae Bloom (HAB) cell counts of *Alexandrium Catenella* (red tide) with toxin concentration (Saxitoxin hydrochloride) in Blue Mussels (*Mytilus edulis*). This project focuses on data collected off of the coast of Maine.

From cell counts, a model of shellfish toxicity is projected over time. Uptake and Depuration constants (alpha and gamma) are determined from correlating historical measured toxicities with modeled toxicity over a range of constants. 

The tools presented herein are designed to present the collected data and calculated toxicity models, and to determine good uptake and depuration constants.

## Installation
Create a conda environment with the following command. Run programs from this git under this shelltox environment.
`conda create -n shelltox python=3.7 bokeh=1.0 pandas=0.23 requests`

## Bokeh App
bokeh_app is a graphical modeling tool for viewing and interacting with data and models.
`restart_bokeh.sh` is a script for restarting web-application service. If installed on a server other that shelltox.whoi.edu it would have to be edited to work properly again.

[img]

## AlphaGamma
alphagamma is a CLI tool for calculating best uptake and depuration constants for a series of years and data-locations. It also calculates estimated best overall constants for a given input series.

Basic usage: `python alphagamma.py alphagamma.ini --src ../data`

Here, `alphagamma.ini` is a configuration file containing year-location triplets with a strong corellation between HAB and measured shellfish toxicity.

