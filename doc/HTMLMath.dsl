<!-- DSSSL Stylesheet fragment HTMLMath.dsl -->

;; Option(s) to pass in relation to LaTeX document class
(define *latexopt* "12pt")

;; LaTeX packages to load?
(define PACKAGES #t)
(define *usepackage* "mathtime")

;; Density specification for the bitmap
(define *density* "96x96")

(root
    (make sequence
      (process-children)
      (process-math)))

(define (write-eqn nd)
  (let ((texmath (select-elements (children (current-node)) 
				  (normalize "alt")))
	(graphic (select-elements (children (current-node)) 
				  (normalize "graphic"))))
    (make element gi: "texequation"
	  attributes: 
	  (list 
	   (list "fileref" (attribute-string (normalize "fileref") graphic)))
	  (literal (data texmath)))))

(mode htmlmath
  (default
    (let ((infeqns (select-elements (descendants (current-node))
				    (normalize "informalequation")))
	  (eqns (select-elements (descendants (current-node))
				 (normalize "equation")))
	  (inleqns (select-elements (descendants (current-node))
				    (normalize "inlineequation"))))
      (with-mode htmlmath
	(process-node-list 
	 (node-list infeqns eqns inleqns)))))

  (element equation (write-eqn (current-node)))
  (element informalequation (write-eqn (current-node)))
  (element inlineequation (write-eqn (current-node))))

(define (process-math)
  (make entity
    system-id: "equation-list.sgml"
    (make element gi: "equation-set"
	  (make element gi: "latexopt"
		(literal *latexopt*))
	  (make element gi: "density"
		(literal *density*))
	  (if PACKAGES  
	      (make element gi: "usepackage"
		    (literal *usepackage*)) (empty-sosofo))
	  (with-mode htmlmath (process-children)))))

