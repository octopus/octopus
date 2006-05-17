#! /usr/bin/python
# -*- coding: utf-8 mode: python -*-
# var2xml.py - A converter to XML for the octopus embeded variabe definitions
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
# $Id: var2xml.py.in 2075 2006-05-14 11:07:56Z athimm $

"""Convert octopus variable definition format into XML.

Takes the variable definitions from stdin and outputs an XML
representation on stdout.

"""

import sys
import fileinput
import re
import cElementTree as ET

root=ET.Element("variables")


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


class LineReader:

    """Convert octopus variable definition format into XML.

    This class is feed linewise with the extracted variable
    definitions from the Fortran source and generates an XML
    representation.

    """

    def __init__(self):
        self.reset()
        
    def reset(self):
        self.variable=""
        self.type=""
        self.sections=""
        self.options=None
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
        theroot=root
        for section in re.split('::',self.sections):
            found=False
            for subelement in theroot.getchildren():
                if subelement.get('name') == section:
                    theroot=subelement
                    found=True
                    break
            if not found:
                theroot=ET.SubElement(theroot,'section',{'name' : section})

        variable=ET.SubElement(theroot,"variable",
                       {"name" : self.variable,
                        "type" : self.type})

        desc = ET.XML("<desc>" + self.desc + "</desc>")
        variable.append(desc)
        if self.options:
            variable.append(self.options)
        if self.default:
            default = ET.XML("<default>" + self.default + "</default>")
            variable.append(default)

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
            if not self.options:
                self.options=ET.Element("options")
            self.optionname=m.group(1)
            self.optionvalue=m.group(2)
            self.multiline=""
            self.readoptiondesc=True
            return
        m=re.match(r'^Option "(.*)"', line)
        if m:
            if not self.options:
                self.options=ET.Element("options")
            self.optionname=""
            self.optionvalue=m.group(1)
            self.multiline=""
            self.readoptiondesc=True
            return
        m=re.match(r'^Option (.*)', line)
        if m:
            if not self.options:
                self.options=ET.Element("options")
            self.optionname=""
            self.optionvalue=m.group(1)
            self.multiline=""
            self.readoptiondesc=True
            return
        if re.match(r'^Description', line):
            self.multiline=""
            return
        print "### ERROR while in " + self.variable + ": ignoring " + line

readme = LineReader()

for line in fileinput.input():
#    sys.stdout.write( "#" + line)
    if line and line.strip():
        readme.feedline(line)

#indent(root)
#
#tree = ET.ElementTree(root)
#tree.write(sys.stdout,"utf-8")

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


printelement(sys.stdout,root)

#ET.ElementTree(root).write("variables-raw.xml")
