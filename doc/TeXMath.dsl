<!-- DSSSL Stylesheet fragment TeXMath.dsl (Allin Cottrell) -->

(element (equation graphic) (empty-sosofo))
(element (equation alt)
 (make display-group
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))

(element (informalequation graphic) (empty-sosofo))
(element (informalequation alt)
 (make display-group
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))
   
(element (inlineequation graphic) (empty-sosofo))
(element (inlineequation alt)
 (make sequence
   (literal "BEGINTEXMATH")
   (literal (data (current-node)))
   (literal "ENDTEXMATH")))


