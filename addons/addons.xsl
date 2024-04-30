<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

<xsl:output method="html" encoding="utf-8"/>

<xsl:template match="/">
  <html>
    <head>
      <title>gretl addons</title>
      <link rel="stylesheet" type="text/css"
      href="../style/server.css"/>
    </head>
    <body>
      <p>[<a href="http://gretl.sourceforge.net/">gretl main page</a>]</p>

      <h1>gretl addons: &#8220;official&#8221; function packages</h1>
      <p>These packages should be installed or updated from within gretl.
      Look for the item &#8220;Check for addons&#8221; under the Help menu.
      </p>
      <p>There are also numerous 
      <a href="https://gretl.sourceforge.net/cgi-bin/gretldata.cgi?opt=SHOW_FUNCS">contributed
       function packages</a> for gretl.</p>
      <div class="gfn-info">
	<xsl:apply-templates/>
      </div>
    </body>
  </html>
</xsl:template>

<xsl:template match="gretl-addon">
  <xsl:call-template name="nl"/>
  <h2><xsl:value-of select="@name"/>: <xsl:value-of select="description"/></h2>
  
  <p>Written by <xsl:value-of select="author"/><br/>
  For gretl <xsl:value-of select="version"/><br/>
  Last revised <xsl:value-of select="date"/><br/>
  JEL code(s): <xsl:value-of select="tags"/><br/>

  Documentation:
  <xsl:element name="a">
    <xsl:attribute name="href">
      <xsl:text>https://sourceforge.net/projects/gretl/files/addons/doc/</xsl:text>
      <xsl:value-of select="helpfile"/>
    </xsl:attribute>
    <xsl:value-of select="helpfile"/>
  </xsl:element>
  </p>
</xsl:template>    

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
