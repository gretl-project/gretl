<!-- DSSSL Stylesheet fragment TeXMath.dsl -->

(element (informalequation texmath)
 (make paragraph
   quadding: 'center
   min-leading: 2pt
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))

(element (informalequation graphic)
  (empty-sosofo))
  

