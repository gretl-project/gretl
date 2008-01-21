<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >

  <!-- Stylesheet for intermediate XML output, for further 
       transformation to TeX -->

<xsl:param name="hlp">cli</xsl:param>
<xsl:param name="standalone">true</xsl:param>
<xsl:param name="lang" select="'en'"/>

<xsl:output method="xml" omit-xml-declaration="yes" indent="yes" encoding="iso-8859-1"/>

<xsl:variable name="intl"
  select="document('hlp_l10n.xml')/internationalization"/>

<xsl:template name="gettext">
  <xsl:param name="key"/>
  <xsl:variable name="itext"
    select="$intl/localization[@language=$lang]/gentext[@key=$key]/@text"/>
  <xsl:choose>
    <xsl:when test="$itext">
      <xsl:value-of select="$itext"/>
    </xsl:when>
    <xsl:otherwise>
      <xsl:message terminate="yes">
        <xsl:text>** Error: no phrase with key = '</xsl:text>
        <xsl:value-of select="$key"/>
        <xsl:text>' found for lang '</xsl:text>
        <xsl:value-of select="$lang"/>
        <xsl:text>'.</xsl:text>
      </xsl:message>
    </xsl:otherwise>
  </xsl:choose>  
</xsl:template>

<xsl:template match="commandlist"> 
 <xsl:apply-templates/> 
</xsl:template>

<xsl:template match="command">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:text>&#xa;</xsl:text>
    <sect2 id="cmd-{@name}" xreflabel="{@name}">
      <title><xsl:value-of select="@name"/></title>
      <xsl:call-template name="nl"/>
      <xsl:apply-templates/>
      <xsl:call-template name="nl"/>
    </sect2>
    <xsl:call-template name="nl"/>    
  </xsl:if>
</xsl:template>

<xsl:template match="description">
  <xsl:if test="not(@context) or @context=$hlp">  
    <xsl:apply-templates/>
  </xsl:if>
</xsl:template>

<xsl:template match="usage">
  <xsl:text>&#xa;</xsl:text>
  <informaltable role="cmd" frame="none">
    <tgroup cols="2"><colspec colnum="1" colwidth="82pt"/>
    <tbody>
      <xsl:call-template name="nl"/>
      <xsl:apply-templates/>
    </tbody>
  </tgroup>
</informaltable> 
<xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="arguments">
 <row>
   <xsl:choose>
     <xsl:when test="count(argument) > 1">
       <entry>
         <para>
           <xsl:call-template name="gettext">
             <xsl:with-param name="key" select="'args'"/>
           </xsl:call-template>
         </para>
       </entry>
     </xsl:when>
     <xsl:otherwise>
       <entry>
         <para>
           <xsl:call-template name="gettext">
             <xsl:with-param name="key" select="'arg'"/>
           </xsl:call-template>
         </para>
       </entry>
     </xsl:otherwise>
   </xsl:choose>
   <entry><para><xsl:apply-templates/></para></entry>
 </row> 
 <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="argblock">
  <xsl:if test="(@optional)"><literal>[ </literal></xsl:if>
  <xsl:if test="(@separated)"><literal>; </literal></xsl:if>
  <xsl:apply-templates/>
  <xsl:if test="(@optional)"><literal>] </literal></xsl:if>
</xsl:template>

<xsl:template match="argument">
  <xsl:if test="(@optional)"><literal>[ </literal></xsl:if>
  <xsl:if test="(@separated)"><literal>; </literal></xsl:if>
  <xsl:if test="(@alternate)">
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'or'"/>
    </xsl:call-template>
  </xsl:if>
  <xsl:if test="@flag">
    <literal><xsl:value-of select="@flag"/></literal>
  </xsl:if> 
  <replaceable><xsl:apply-templates/></replaceable>
  <xsl:text> </xsl:text>
  <xsl:if test="(@optional)"><literal>] </literal></xsl:if>
</xsl:template>

<xsl:template match="argpunct">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="altform">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <entry>
          <para>
            <xsl:call-template name="gettext">
              <xsl:with-param name="key" select="'altforms'"/>
            </xsl:call-template>
          </para>
        </entry>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><para><xsl:apply-templates/></para></entry>
  </row>
  <xsl:text>&#xa;</xsl:text>
</xsl:template>

<xsl:template match="option">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../option) > 1">
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'opts'"/>
                </xsl:call-template>
              </para> 
            </entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'opt'"/>
                </xsl:call-template>
              </para> 
            </entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><para><xsl:apply-templates/></para></entry>
  </row>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="example">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../example) > 1">
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'examples'"/>
                </xsl:call-template>
              </para>
            </entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'example'"/>
                </xsl:call-template>
              </para>
            </entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry><para><literal><xsl:apply-templates/></literal></para></entry>
  </row>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demos">
  <row>
    <xsl:choose>
      <xsl:when test="position() = 1">
        <xsl:choose>
          <xsl:when test="count(../example) > 1">
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'examples'"/>
                </xsl:call-template>
              </para>
            </entry>
          </xsl:when>
          <xsl:otherwise>
            <entry>
              <para>
                <xsl:call-template name="gettext">
                  <xsl:with-param name="key" select="'example'"/>
                </xsl:call-template>
              </para>
            </entry>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:otherwise>
        <entry/>
      </xsl:otherwise>
    </xsl:choose>
    <entry>
      <para>
        <xsl:if test="count(../example) > 1">
          <xsl:call-template name="gettext">
            <xsl:with-param name="key" select="'Seealso'"/>
          </xsl:call-template>
        </xsl:if>
        <xsl:apply-templates/>
      </para>
    </entry>
  </row>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="demo">
  <xsl:if test="position() > 1">
    <xsl:text>, </xsl:text>
  </xsl:if>
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="flag">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="effect">
  <xsl:text> (</xsl:text><xsl:apply-templates/><xsl:text>)</xsl:text>
</xsl:template>

<xsl:template match="repl">
  <replaceable><xsl:apply-templates/></replaceable>
</xsl:template>

<xsl:template match="cmd">
  <command><xsl:apply-templates/></command>
</xsl:template>

<xsl:template match="program">
  <application><xsl:apply-templates/></application>
</xsl:template>

<xsl:template match="lit">
  <literal><xsl:apply-templates/></literal>
</xsl:template>

<xsl:template match="mathvar">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>

<xsl:template match="emphasis">
  <emphasis><xsl:apply-templates/></emphasis>
</xsl:template>

<xsl:template match="super">
  <superscript><xsl:apply-templates/></superscript>
</xsl:template>

<xsl:template match="sub">
  <subscript><xsl:apply-templates/></subscript>
</xsl:template>

<xsl:template match="book">
  <citetitle><xsl:apply-templates/></citetitle>
</xsl:template>

<xsl:template match="quote">
  <quote><xsl:apply-templates/></quote>
</xsl:template>

<xsl:template match="footnote">
  <footnote><xsl:apply-templates/></footnote>
</xsl:template>

<xsl:template match="filename">
  <filename><xsl:apply-templates/></filename>
</xsl:template>

<xsl:template match="func">
  <function><xsl:apply-templates/></function>
</xsl:template>

<xsl:template match="equation">
 <xsl:choose>
  <xsl:when test="(@status='display')">
   <informalequation>
    <alt role="tex"><xsl:value-of select="@tex"/></alt>
    <graphic fileref="figures/{@graphic}"/>
   </informalequation>
  </xsl:when>
  <xsl:otherwise>
   <inlineequation>
    <alt role="tex"><xsl:value-of select="@tex"/></alt> 
    <inlinemediaobject>
     <imageobject>
      <imagedata align="center" fileref="figures/{@graphic}"/>
     </imageobject>
    </inlinemediaobject>
   </inlineequation>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template match="ilist">
  <xsl:if test="not(@context) or @context=$hlp">
    <itemizedlist><xsl:apply-templates/></itemizedlist>
    <xsl:call-template name="nl"/>
  </xsl:if>
</xsl:template>

<xsl:template match="nlist">
  <xsl:if test="not(@context) or @context=$hlp">
    <orderedlist><xsl:apply-templates/></orderedlist>
    <xsl:call-template name="nl"/>
  </xsl:if>  
</xsl:template>

<xsl:template match="li">
  <listitem><xsl:apply-templates/></listitem>
  <xsl:call-template name="nl"/>
</xsl:template>

<xsl:template match="para">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:call-template name="nl"/>
    <para>
      <xsl:apply-templates/>
    </para>
    <xsl:call-template name="nl"/>
  </xsl:if>   
</xsl:template>

<xsl:template match="pre">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:call-template name="nl"/>
    <para role="preformatted">
      <xsl:apply-templates/>
    </para>
    <xsl:call-template name="nl"/>
  </xsl:if>   
</xsl:template>

<xsl:template match="code">
  <xsl:if test="not(@context) or @context=$hlp">
    <programlisting>
      <xsl:apply-templates/>
    </programlisting>
   </xsl:if>   
</xsl:template>

<xsl:template match="cmdref">
  <xref linkend="cmd-{@targ}"/>
</xsl:template>

<xsl:template match="guideref">
  <xsl:choose>
    <xsl:when test="$standalone='true'">
      <xsl:text>\GUG{}</xsl:text>
    </xsl:when>
    <xsl:otherwise>
      <xref linkend="{@targ}"/>
    </xsl:otherwise>
  </xsl:choose>
</xsl:template>

<xsl:template match="menu-path">
  <para>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'menupath'"/>
    </xsl:call-template>
  <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="other-access">
  <para>
    <xsl:call-template name="gettext">
      <xsl:with-param name="key" select="'otheraccess'"/>
    </xsl:call-template>
  <xsl:apply-templates/></para>
</xsl:template>

<xsl:template match="table">
  <xsl:if test="not(@context) or @context=$hlp">
    <xsl:choose>
      <xsl:when test="@id">
        <table id="{@id}" frame="none">
          <title><xsl:value-of select="@title"/></title>
          <tgroup cols="2" style="{@style}">
            <xsl:if test="@lwidth and @rwidth">
              <colspec colwidth="{@lwidth}"/>
              <colspec colwidth="{@rwidth}"/>
            </xsl:if>
            <thead>
              <row>
                <entry><xsl:value-of select="@lhead"/></entry>
                <entry><xsl:value-of select="@rhead"/></entry>
              </row>
            </thead>
            <tbody>
              <xsl:apply-templates/>
            </tbody>
          </tgroup>
        </table>
      </xsl:when>
      <xsl:otherwise>
        <informaltable role="cmd" frame="none">
          <tgroup cols="2" style="{@style}">
            <xsl:if test="@lwidth and @rwidth">
              <colspec colwidth="{@lwidth}"/>
              <colspec colwidth="{@rwidth}"/>
            </xsl:if>
            <tbody>
              <xsl:apply-templates/>
            </tbody>
          </tgroup>
        </informaltable>
      </xsl:otherwise>
    </xsl:choose>    
  </xsl:if>
</xsl:template>

<xsl:template match="row">
  <row><xsl:apply-templates/></row>
</xsl:template>

<xsl:template match="cell">
  <entry><para><xsl:apply-templates/></para></entry>
</xsl:template>

<xsl:template name="nl">
  <xsl:text>&#10;</xsl:text>  
</xsl:template>

</xsl:stylesheet>
