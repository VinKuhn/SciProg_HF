# SciProg_HF
Hartree Fock program for calculations with the restricted Hartree-Fock model coded in Fortran
## Table of Contents
* [Introduction](#introduction)
* [General info](#general-info)
* [Input](#input)
* [Technologies](#technologies)
* [Setup](#setup)

## Introduction
This project was based on a provided skeleton HF program provided by Luuk Visscher, that used only a limited set of base functions for a preset Be-He molecule. The incentive was to extend the program to accept input for any molecule and also implement the iterative process to create the Fock matrix. 

## General info
Project is created with:
* Fortran90
* gen1int library
* SciProg_HF by Luuk Visscher

## Input
The input file has to be named molecule.xyz and provided within the projects folder, also x, y and z coordinates of the atoms have to be provided in atomic units. Formatting of the input is as follows:



![template](https://user-images.githubusercontent.com/101809431/161159992-ad7848bb-eec2-4f2d-8679-a9a6ce1ed844.png)


Following that scheme the input for a He-He molecule with a distance of 2 a.u. would look as follows:


![input](https://user-images.githubusercontent.com/101809431/161159744-d6ce5bf9-ffc9-4b0d-9100-a5805e7b456a.png)

!Please note, that reasonable results will only be achhieved for molecules with an even number of electrons (restricted HF)


