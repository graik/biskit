Biskit
======

Main Classes
------------

The most important classes which can be directly imported from the `biskit`
namespace.

.. currentmodule:: biskit
                   
.. autosummary::
   :nosignatures:
   :toctree: biskit
   :template: autosummary/class_with_details.rst       

   PDBModel
   ProfileCollection   
   PDBCleaner
   Trajectory
   EnsembleTraj
   Residue
   AmberResidueType
   AmberPrepParser
   LogFile
   StdLog
   ErrLog

Interaction with external programs
----------------------------------

.. autosummary::
   :nosignatures:
   :toctree: biskit
   :template: autosummary/class_with_details.rst       

    Executor
    ExeConfig
    ExeConfigCache
    
    AmberCrdEntropist
    AmberCrdParser
    AmberLeap
    AmberParmBuilder
    AmberRstParser
    AmberResidueLibrary
    
    TMAlign
    
Supporting classes and modules
------------------------------

.. autosummary::
   :nosignatures:
   :toctree: biskit
   :template: autosummary/class_with_details.rst
   
   BioUnit
   LocalPath
   SettingsManager
   SettingsParser


              
Helper Modules
--------------

Modules with utility and helper methods. 

.. currentmodule:: biskit

.. autosummary::
   :toctree: biskit
   :template: autosummary/module_with_details.rst
                        
   gnuplot
   hist
   tools
   molUtils
   mathUtils
   match2seq           
   rmsFit
   
Errors
-------

Errors that are raised by biskit classes.

.. currentmodule:: biskit

.. autosummary::
   :nosignatures:
   :toctree: biskit
   :template: autosummary/class_with_details.rst       
   
   EHandler
   BiskitError
   LocalPathError
   PDBError
   ProfileError
