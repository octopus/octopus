#! @PYTHON@
# -*- coding: utf-8 mode: python -*-
# XML.py - Manage XML input and output for octopus variables.
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

"""Manage XML input and output for octopus variables.

This modules takes care of writing out variable definitions and values
into files and reading them back."""

import base
import re
import StringIO
import re

import cElementTree as ET

import base

def fixxmlstring(string):
    """Return XML-sanitized string.
    
    This method replaces all bugus (from XML's POV) entities. It
    later restores 'known' tags. Finally it also adds newlines
    after <br />.
        
    """

    string = string.replace('&', '&amp;')
    string = string.replace('"', '&quot;')
    #string = string.replace('<=', '&le;')
    #string = string.replace('>=', '&ge;')
    string = string.replace('<', '&lt;')
    string = string.replace('>', '&gt;')
    for tag in ['tt', 'math' ,'a', 'br', 'i', 'li', 'ul']:
        string = re.sub(r'&lt;' + tag + r'&gt;', '<' + tag + r'>', string)
        string = re.sub(r'&lt;' + tag + r'( .*?)&gt;', '<' + tag + r'\1>', string)
        string = string.replace('<br>', '<br />')
        string = string.replace('&lt;/' + tag + '&gt;', '</' + tag + '>')
    return string

def var2element(variable):
    element=ET.Element("variable",
                          {"name" : variable.name,
                           "type" : variable.type})
    element.append(ET.XML("<desc>" + variable.description + "</desc>"))
    if "options" in variable.__dict__:
        optelement=ET.Element("options")
        for option in variable.options:
            optelement.append(option)
        element.append(optelement)
    if "default" in variable.__dict__:
        element.append(ET.XML("<default>" + variable.default + "</default>"))
    return element


def element2var(element):
    name=element.get("name")
    type=element.get("type")
    description="no description available"
    for child in element.getchildren():
        if child.tag=="desc" and ((child.text and child.text.strip())
                                  or child.getchildren() != []):
                strio=StringIO.StringIO()
                #ET.ElementTree(child).write(strio)
                printelement(strio,child)
                m=re.match(r'<desc>\s*(.*?)\s*</desc>', strio.getvalue(), re.S | re.M)
                description=m.group(1)
    return base.variable_by_type(name,type,description)

def MakeElement(variabletree, root=ET.Element("variables")):

    for variable in variabletree.variables:
        root.append(var2element(variable))
    for section in variabletree.subsections:
        MakeElement(section,
                    ET.SubElement(root, "section", {"name" : section.name}))
    return root

class LineReader(object):

    """Convert octopus variable definition format into XML.

    This class is feed linewise with the extracted variable
    definitions from the Fortran source and generates an XML
    representation.

    """

    def __init__(self,variabletree):
        self.reset()
        self.variabletree=variabletree
        
    def reset(self):
        self.variable=""
        self.type=""
        self.sections=""
        self.options=[]
        self.optionname=""
        self.optionvalue=""
        self.default=""
        self.desc=""
        self.multiline=None
        self.readoptiondesc=False

    def writeoutvariable(self):
        """Create a variable element.

        This method collates all parsed information about the variable
        into an element and adds the element within the main tree. The
        'section' information is used to determine where in the tree
        the element will be appended to. In case the section tree part
        is missing it is created on the fly.

        """

        # Find/Create the (sub)section
        theroot=self.variabletree
        for section in re.split('::',self.sections):
            found=False
            for subelement in theroot.get_sections():
                if subelement.name == section:
                    theroot=subelement
                    found=True
                    break
            if not found:
                theroot=theroot.add_section(base.VariableTree(section))

        if self.options != []:
            self.type="option"

        variable=base.variable_by_type(self.variable,self.type,self.desc)
        theroot.add_variable(variable)

        if self.options != []:
            variable.options=self.options

        if self.default:
            default = ET.XML("<default>" + self.default + "</default>")
            variable.default=self.default

        self.reset()

    def feedline(self,line):
        """Parse line for XML conversion.

        This method takes a complete line as input and dependent on
        the mode it is in (normal or multiline) it will either
        continue constucting a multiline input (variable or options'
        descriptions) or call parseline. Special handling happens for
        the 'End' tag.
        
        """

        if self.multiline != None:
            if re.match(r'^ ', line):
                self.multiline += ' ' + line.strip()
                return
            else:
                self.multiline = fixxmlstring(self.multiline.strip())
                if(self.readoptiondesc):
                    if self.optionname=="":
                        option=ET.XML('<option value="' + self.optionvalue + '">'
                                      + self.multiline + '</option>')
                    else:
                        option=ET.XML('<option name="' + self.optionname + '" value="' + self.optionvalue + '">'
                                      + self.multiline + '</option>')
                    self.options.append(option)
                else:
                    self.desc=self.multiline
                self.multiline=None
                self.readoptiondesc=False
        if re.match(r'^End.*', line):
            self.writeoutvariable()
            return
        else:
            self.parseline(line)

    def parseline(self,line):
        """Parse non-multiline and no-'End' constucts.

        This method takes a line assuming that it is not a multiline
        and parses it. It doesn't expect 'End' tags (these are taken
        care of feedline). If an unknown tag is found it gives an
        error, but doesn't stop processing.

        """
        m=re.match(r'^Variable (.*)', line)
        if m:
            self.variable=m.group(1)
            return
        m=re.match(r'^Type (.*)', line)
        if m:
            self.type=m.group(1)
            return
        m=re.match(r'^Section (.*)', line)
        if m:
            self.sections=m.group(1)
            return
        m=re.match(r'^Default (.*)', line)
        if m:
            self.default=m.group(1)
            return
        m=re.match(r'^Option (.*) (.*)', line)
        if m:
            self.optionname=m.group(1)
            self.optionvalue=m.group(2)
            self.multiline=""
            self.readoptiondesc=True
            return
        m=re.match(r'^Option "(.*)"', line)
        if m:
            self.optionname=""
            self.optionvalue=m.group(1)
            self.multiline=""
            self.readoptiondesc=True
            return
        m=re.match(r'^Option (.*)', line)
        if m:
            self.optionname=""
            self.optionvalue=m.group(1)
            self.multiline=""
            self.readoptiondesc=True
            return
        if re.match(r'^Description', line):
            self.multiline=""
            return
        print "### ERROR while in " + self.variable + ": ignoring " + line

def isblockelement(element):
    return element.tag in ('variables', 'section', 'variable', 'desc', 'options', 'option', 'ul', 'li', 'default')

def iscompactblockelement(element):
    return element.tag in ('li', 'default')

def printelement(out,element,level=0,blockout=True):
    hastext = element.text and element.text.strip()
    hastail = element.tail and element.tail.strip()
    children = element.getchildren()
    haschildren = (children != [])
    blockout = blockout and isblockelement(element)
    blockin = blockout and not iscompactblockelement(element)
    tags=""
    for attr in element.items():
        tags += ' ' + attr[0] + '="' + attr[1] + '"'
    # element with content?
    if hastext or haschildren:
        # start tag
        out.write('<' + element.tag + tags + '>')
        # text
        if hastext:
            if blockin:
                out.write('\n' + (level+1)*'  ')
            out.write(fixxmlstring(element.text))
        # children
        if haschildren:
            for child in children:
                if blockin and isblockelement(child):
                    out.write('\n' + (level+1)*'  ')
                printelement(out,child,level+1,blockin)
        # end tag
        if blockin:
            out.write('\n' + level*'  ')
        out.write('</' + element.tag + '>')
#        if block:
#            out.write('\n')
    else: # has no content
        out.write('<' + element.tag + tags + ' />')
        if element.tag == 'br':
            out.write('\n' + (level-1)*'  ')
    if hastail:
        out.write(fixxmlstring(element.tail))

