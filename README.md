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
The input file has to be named molecule.xyz and provided within the projects folder. 
The x, y and z coordinates have to be provided in atomic units
Formatting of the input is as follows:

n_atoms        molecular_charge
nuclear_charge1      x             y            z
nuclear_charge2      x             y            z
      .              .             .            .
      .              .             .            .
Following that scheme the input for a He-He molecule with a distance of 2 a.u. would look as follows:

*2 0
*2 0 0 0
*2 2 0 0

!Please note, that reasonable results will only be achhieved for molecules with an even number of electrons (restricted HF)


