<?xml version='1.0'?>

<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

<xsl:template match="/">
  <html>
    <xsl:apply-templates/>
  </html>
</xsl:template>

<xsl:template match="config_compset">
  <head>
    <title>Configuration Component Sets</title>
  </head>
  <body>
    <h2>Configuration Component Sets</h2>

    <table BORDER="1" CELLPADDING="10">
      <th>Name (Shortname)</th>
      <th>Description, Including components, CCSM variables</th>
      <tr>
      <th colspan="3">A (All Data Models)</th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'A_')]"/>
      <tr>
      <th colspan="3">B (All Active Models)</th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'B_')]"/>
      <tr>
      <th colspan="3">C (Standalone POP)</th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'C_')]"/>
      <tr>
      <th colspan="3">D (Active sea-ice and ocean, data atmosphere and stub land)</th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'D_')]"/>
      <tr>
      <th colspan="3">E (Active land and atmosphere with slab ocean) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'E_')]"/>
      <tr>
      <th colspan="3">F (Active land and atmosphere with data ocean) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'F_')]"/>
      <tr>
      <th colspan="3">G (Active sea-ice and ocean, data atmosphere and land) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'G_')]"/>
      <tr>
      <th colspan="3">H (Active sea-ice and ocean, data atmosphere and stub land) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'H_')]"/>
      <tr>
      <th colspan="3">I (Standalone CLM) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'I_')]"/>
      <tr>
      <th colspan="3">S (All stub models with dead atmosphere model) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'S_')]"/>
      <tr>
      <th colspan="3">X (All dead models) </th>
      </tr>
      <xsl:apply-templates select="compset[contains(@NAME,'X')]"/>
    </table>

  </body>
</xsl:template>

<xsl:template match="compset">
  <tr>
    <td><font color="#ff0000"><xsl:value-of select="@NAME"/></font>
    (<xsl:value-of select="@SHORTNAME"/>)
        <xsl:if test="string-length(@GRID_MATCH)>0">
        [Grid=<xsl:value-of select="@GRID_MATCH"/>]
        </xsl:if>
    </td>
    <td><xsl:value-of select="@DESC"/>
    <xsl:if test="string-length(@COMP_ATM)>0">
    <p>
    (atm=<xsl:value-of select="@COMP_ATM"/>
     lnd=<xsl:value-of select="@COMP_LND"/>
     glc=<xsl:value-of select="@COMP_GLC"/>
     ice=<xsl:value-of select="@COMP_ICE"/>
     ocn=<xsl:value-of select="@COMP_OCN"/>)
    </p>
    </xsl:if>
        <xsl:if test="string-length(@CCSM_CO2_PPMV)>0">
        CCSM_CO2_PPMV=<xsl:value-of select="@CCSM_CO2_PPMV"/>
        </xsl:if>
        <xsl:if test="string-length(@CCSM_BGC)>0">
        CCSM_BGC=<xsl:value-of select="@CCSM_BGC"/>
        </xsl:if>
        <xsl:if test="string-length(@RUN_REFCASE)>0">
        ref_case=<xsl:value-of select="@RUN_REFCASE"/>
        </xsl:if>
        <xsl:if test="string-length(@RUN_REFCASE)>0">
        ref_date=<xsl:value-of select="@RUN_REFDATE"/>
        </xsl:if>

    </td>
  </tr>
</xsl:template>


</xsl:stylesheet>
