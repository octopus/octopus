#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# base.by - Basic definitions for octopus variables.
# Copyright (C) 2006 A. Thimm
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
# 
# $Id$

"""Basic definitions for octopus variables.

This module contains class definitions and basic operation for octopus
variables."""

# FIXME: merge variable sets
# FIXME: differentiate between definitions and values?

class OctopusVariable(object):

    def __init__(self, name, type, description):
        self.name = name
        self.type = type
        self.description = description

class OctopusVariableFloat(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "float", description)

class OctopusVariableInteger(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "integer", description)

class OctopusVariableOption(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "option", description)

class OctopusVariableString(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "string", description)

class OctopusVariableLogical(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "logical", description)

class OctopusVariableFlag(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "flag", description)

class OctopusVariableBlock(OctopusVariable):

    def __init__(self, name, description):
        OctopusVariable.__init__(self, name, "block", description)

class Section(object):

    def __init__(self, name, description):
        self.subsections=[]
        self.variables=[]
        self.name=name
        self.description=description

def variable_by_type(name,type,description):

    var=None
    if type=="float": var=OctopusVariableFloat(name, description)
    if type=="integer": var=OctopusVariableInteger(name, description)
    if type=="option": var=OctopusVariableOption(name, description)
    if type=="string": var=OctopusVariableString(name, description)
    if type=="logical": var=OctopusVariableLogical(name, description)
    if type=="flag": var=OctopusVariableFlag(name, description)
    if type=="block": var=OctopusVariableBlock(name, description)
    return var

class VariableTree(object):

    def __init__(self, name="", description=""):
        self.name=name
        self.description=description
        self.subsections=[]
        self.variables=[]

    def add_section(self, section):
        self.subsections.append(section)
        return section

    def get_sections(self):
        return self.subsections

    def add_variable(self, variable):
        self.variables.append(variable)

    def get_variables(self):
        return self.variables

    def print2(self,out):
        for variable in self.variables:
            out.write(variable.name)
        for section in self.subsections:
            out.write(section.name)
            section.print2(out)
            
