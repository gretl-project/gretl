<!-- DSSSL Stylesheet fragment HTMLMath.dsl -->

;; User-configurable options

(define $latexopt$ 
  ;; Option(s) to pass in relation to LaTeX article document class
  ;; Set to empty string "" for no option
  "12pt")

(define $usepackage$
  ;; LaTeX packages to load?
  ;; Set to empty string "" for no packages
  ;; Can use, e.g., "mathtime" to use MathTime fonts
  ;; "mathtime")  
  "mathtime")

(define $density$
  ;; Density specification for equation bitmaps
  ;; This will be passed to the "convert" program
  "96x96")

;; End of user-configurable options  

;; stuff to avoid leading newlines messing up the TeX when using
;; xml docbook.

(define %newline% "
")

(define (count-newlines str)
  (let loop ((count 0))
    (if (string=? (substring str count (+ count 1)) %newline%) 
        (loop (+ count 1))
        count
    )
  ))

;; Returns string, with any leading newlines trimmed off
(define (trim-leading string)
  (substring string
             (count-newlines string)
             (string-length string)))

(root
    (make sequence
      (process-children)
      (process-math)))

;; How to write out an equation into the equation listing file
(define (write-eqn nd)
  (let ((texmath (select-elements (children (current-node)) 
				  (normalize "alt")))
	(graphic (select-elements (children (current-node)) 
				  (normalize "graphic"))))
    (make element gi: "texequation"
	  attributes: 
	  (list 
	   (list "fileref" (attribute-string (normalize "fileref") graphic)))
	  (literal (trim-leading (data texmath))))))

;; Special processing mode to extract equations
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

;; Write equation info to equation-list.sgml
(define (process-math)
  (make entity
    system-id: "equation-list.sgml"
    (make element gi: "equation-set"
	  attributes: (list
		       (list "latexopt" $latexopt$)
		       (list "density" $density$)
		       (list "usepackage" $usepackage$))
	  (with-mode htmlmath (process-children)))))
