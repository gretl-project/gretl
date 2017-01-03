# -*- coding: utf-8 -*-

import re
import libxml2
import os
import sys
try:
    # Hashlib is new in Python 2.5
    from hashlib import md5 as md5_new
except ImportError:
    from md5 import new as md5_new

from basic import basicXmlMode

class gretldocXmlMode(basicXmlMode):
    """Class for special handling of gretl document types.

    It sets lang attribute on article elements, and adds translators
    to articleinfo/copyright."""
    def __init__(self):
        self.lists = ['ilist', 'nlist']
        self.objects = ['seelist', 'equation', 'code', 'math']

    def getIgnoredTags(self):
        "Returns array of tags to be ignored."
        return self.objects + self.lists

    def getFinalTags(self):
        "Returns array of tags to be considered 'final'."
        return ['para', 'fnarg'] + self.objects + self.lists

    def getSpacePreserveTags(self):
        "Returns array of tags in which spaces are to be preserved."
        return ['tabular']

    def preProcessXml(self, doc, msg):
        "Preprocess a document and perhaps adds some messages."
        pass

    def postProcessXmlTranslation(self, doc, language, translators):
        """Sets a language and translators in "doc" tree.

        "translators" is a string consisted of translator credits.
        "language" is a simple string.
        "doc" is a libxml2.xmlDoc instance."""
        pass

    def getStringForTranslators(self):
        """Returns None or a string to be added to PO files.

        Common example is 'translator-credits'."""
        return None

    def getCommentForTranslators(self):
        """Returns a comment to be added next to string for crediting translators.

        It should explain the format of the string provided by getStringForTranslators()."""
        return None

# Perform some tests when ran standalone
if __name__ == '__main__':
    test = gretlXmlMode()
    print "Ignored tags       : " + repr(test.getIgnoredTags())
    print "Final tags         : " + repr(test.getFinalTags())
    print "Space-preserve tags: " + repr(test.getSpacePreserveTags())

