<!-- DSSSL Stylesheet fragment mathml.dsl
     (included as an entity into mathmlx.dsl)
  -->

<!--docbookhack set -> mset -popel-->
<!--revised by Jan Scheffczyk (herta@mail.studfb.unibw-muenchen.de)-->

;;;
;;;    David Carlisle
;;;    davidc@nag.co.uk
;;;
;;;    Copyright 1998 Nag Ltd, The OpenMath Consortium. Esprit Project 24.969.
;;;



;;;     COLOUR
;;;
;;; Just RGB colour supported currentlly

(define (rgb-color r g b)
           (color (color-space 
	      "ISO/IEC 10179:1996//Color-Space Family::Device RGB") r g b))



;;;   CONTENT MathML
;;;
;;;   Content MathML is mainly implemented directly with element
;;;   declarations, and process-children. there is not too much
;;;   need for node list processing.
;;;
;;;   Currently many attributes for font changes and spacing are silently
;;;   ignored.
;;;
;;;   The mo element goes to some trouble to get its attributes as specified
;;;   in the MathML recomendation, but currently doesn't do much with them.


;;; mrow
;;; should check attributes (this comment applies to most elements
;;; but won't be repeated

(element mml:mrow
  (process-children-trim))


;;; mi
;;; Math Identifier Defaults to italic.
;;; Ought to switch between math italic and text italic
;;; for multi letter identifiers (or just in tex backend?)

(element mml:mi
  (make math-sequence
    font-posture: 
    (let ((fnt 
	   (if (attribute-string "fontstyle")
	      (attribute-string "fontstyle")
	      "italic")))
      (if (equal? "style-normal" fnt)
	  'upright
	  (if (equal? "italic" fnt)
	      'italic 
	      #f)))
    (process-children-trim)))


;;; mn
;;; Same for numbers.

(element mml:mn
  (make math-sequence
  font-posture: 
  (let ((fnt 
       (if(attribute-string "fontstyle")
        (attribute-string "fontstyle")
           "style-normal")))
    (if (equal? "style-normal" fnt)
        'upright
        (if (equal? "italic" fnt)
          'italic 
            #f)))
        (process-children-trim)))

;;; mtext
;;; Bits of non-math

(element mml:mtext
  (make unmath
   (process-children-trim)))

;;; mspace
;;; Grumble grumble it seems extraordinarily complicated to copy
;;; a length from an attribute on the element to a keyword to a make
;;; function. Also the rtf backend doesn't really support line-field
;;; I couldn't get inline-space characters to work either.

(element mml:mspace
  (make line-field   field-width: 
         (let ((x   (attribute-value "width" (current-node))));;;???
         (measurement-to-length (if (attribute-string "width")
        (attribute-string "width")
           "0pt")))))

;;; ms
;;; Doesn't work right in tex backend: How do you specify mono space font
;;; without specifying what font to use.

(element mml:ms
  (make unmath
    font-posture: 'upright
    font-family-name: "iso-monospace"
      (literal "\"")
      (process-children-trim)
      (literal "\"")))


;;; mfrac
;;; fractions.

(element mml:mfrac
  (make fraction
    (let ((nl (children(current-node))))
     (sosofo-append
      (make math-sequence
        label: 'numerator
         (process-node-list (node-list-first nl)))
      (make math-sequence
        label: 'denominator
         (process-node-list (node-list-rest nl)))))))

;;; msqrt mroot
;;; Radicals

(element mml:msqrt
  (make radical
   (process-children-trim)))


(element mml:mroot
  (make radical
    (let ((nl (children(current-node))))
     (sosofo-append
      (make math-sequence
         (process-node-list (node-list-first nl)))
      (make math-sequence
        label: 'degree
         (process-node-list (node-list-rest nl)))))))


;;; mstyle
;;; Style, what style?

(element mml:mstyle
  (make math-sequence
     (process-children-trim)))

;;; merror
;;; Ignore this, for now

(element mml:merror
  (make math-sequence
    (process-children-trim)))

;;; mpadded
;;; Hmmm
(element mml:mpadded
  (make math-sequence
    (process-children-trim)))


;;; mphantom
;; do it in white: not really the same as invisible
;; but not sure if there is an easy general way to access
;; background colour.

(element mml:mphantom
  (make math-sequence
     color: (rgb-color  1 1 1)
    (process-children-trim)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; mfenced
;;; Doesn't do separators for now.

(element mml:mfenced
  (make fence
      (if  (attribute-string "open")
      (make math-sequence
         label: 'open
          (literal (attribute-string "open"))) 
      (empty-sosofo))
      (if  (attribute-string "close")
        (make math-sequence
         label: 'close
          (literal (attribute-string "close")))
         (empty-sosofo))
    (make math-sequence
    (process-children-trim))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Superscripts and subscripts

;;; msup

(element mml:msup
  (make script
    (let ((nl (children(current-node))))
     (sosofo-append
      (make math-sequence
         (process-node-list (node-list-first nl)))
      (make math-sequence
        label: 'post-sup
         (process-node-list (node-list-rest nl)))))))

;;; msub

(element mml:msub
  (make script
    (let ((nl (children(current-node))))
     (sosofo-append
      (make math-sequence
         (process-node-list (node-list-first nl)))
      (make math-sequence
        label: 'post-sub
         (process-node-list (node-list-rest nl)))))))

;;; msubsup

(element mml:msubsup
  (make script
    (let* ((nl (children(current-node)))
              (nlr (node-list-rest nl)))
     (sosofo-append
      (make math-sequence
            (process-node-list (node-list-first nl)))
      (make math-sequence
         label: 'post-sub
         (process-node-list (node-list-first nlr)))
      (make math-sequence
        label: 'post-sup
         (process-node-list (node-list-rest nlr)))))))

;;; mmultiscripts
;;; In order to get the scripts aligning with each other
;;; they all script an empty element (so ignore the size of the base
;;; I wish I could measure things in DSSSL....

(element mml:mmultiscripts
    (let* ((nl (children(current-node)))
              (base (node-list-first nl))
              (nlr (node-list-rest nl)))
      (process-multi-scripts base nlr #t (empty-sosofo) (empty-sosofo))))

  
;;; while flag is true scoop up the scripts into the fourth argument
;;; then when you see multiscripts switch the flag so then start collecting
;;; in the third argument. Finally when rest is empty, stuff the scripts
;;; around the base.
       
(define (process-multi-scripts base rest flag left right )
  (if (node-list-empty? rest)
      (sosofo-append 
       left
       (make math-sequence 
	 (process-node-list base))
       right)
; else
      (if (equal? "MPRESCRIPTS" (gi (node-list-first rest)))
	  (process-multi-scripts base (node-list-rest rest) #f left right)
	  (if flag
	      (process-multi-scripts
	       base
	       (node-list-rest (node-list-rest rest))
	       flag 
	       left
	       (sosofo-append
		right
		(make script
		  (make math-sequence
		    label: 'post-sub
		    (process-node-list (node-list-first rest)))
		  (make math-sequence
		    label: 'post-sup
		    (process-node-list (node-list-first (node-list-rest rest)))))))
	; else
	      (process-multi-scripts
	       base
	       (node-list-rest (node-list-rest rest))
	       flag
	       (sosofo-append
		left
		(make script
		  (make math-sequence
		    label: 'post-sub
		    (process-node-list (node-list-first rest)))
		  (make math-sequence
		    label: 'post-sup
		    (process-node-list (node-list-first (node-list-rest rest))))))
	       right)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; none/
;;; Not the hardest thing to implement

(element none
  (empty-sosofo))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Under and over bars.
;;; These don't work yet.

(element mml:munder
  (make math-sequence
   (process-children-trim)))

(element mml:mover
  (make math-sequence
   (process-children-trim)))

(element mml:munderover
  (make math-sequence
   (process-children-trim)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Operator dictionary

;;; This is just a small version, while I test various implementations.


;;; First put back the indirection that TeX as but that MathML omitted.
;;; Need the default lengths to refer to variables rather than explicit
;;; lengths, so that you can change them all together.
;;; Apart from the fact that it isn't implemented, there arelots of other
;;; problems with this operator dictionay concept.

;;; (mathml problem) As noted above the recommendation gives explicit lengths.
;;; (mathml problem) As the recommendation gives fixed (as opposed to
;;;                  variable/stretchy) lengths.
;;; (dsssl problem)  Full support of stretch operators would require an
;;;                  extended dsssl flow object (or extended semantics
;;;                  for the stretchy character characteristic) to support
;;;                  `mid' fences.
;;; 
(define %zskip 0em)
(define %smallskip .05555em)
(define %medskip .16666em)
(define %bigskip .16666em)

(define mathml-op-dict
'(
 ("+" 
   ("prefix". ((lspace . %zskip) (rspace . %bigskip)))
   ("infix" . ((lspace . %medskip) (rspace . %medskip))))
 ("*="
   ("infix" .  ((lspace . %bigskip)(rspace . %bigskip))))
 ("("
   ("prefix" . ((fence . #t )(stretchy . #t)(lspace . %zskip)(rspace . %zskip))))
 (")"
   ("postfix" .((fence . #t )(stretchy . #t)(lspace . %zskip)(rspace . %zskip))))

))


;;; Helper function, just trims spaces from strings, actually
;;; zaps all space. needed as mo doesn't use process-children-trim.

(define (string-nospace s)
  (let ((l (string-length s)))
    (let loop ((i 0))
      (if (= i l)
          ""
          (let (( x (string-ref s i)))
          (if (equal? #\space x)
              (loop (+ i 1))
              (string-append (string x)
                (loop (+ i 1)))))))))

;;; mo
;;; This current version tries to implement the defaulting
;;; of the various attributes. Although currently
;;; It doesn't actually use the attributes once it has got the
;;; values.

(element mml:mo
 (make math-sequence

  (let* 
;; First get the name of the operator
;; and the surrounding node list which will be dealt with specially
;; if it is mrow (or more correctly should be mrow-like)
      ((nm (string-append(string-nospace (data (current-node)))))
       (pnt (parent))
;; An explicit form attribute (or #f)
       (form1 (attribute-string "form"))
;; Look up the name in the operator dictionary, which will return
;; another association list, for looking up the form [if the operator
;; is in the dictionary].
       (p (assoc nm mathml-op-dict))
       (form
	(if form1
	    (if p
;; If a form was specified, and the operator is in the dictionary
;; look up the list of defaults. If it isn't in the dictionary
;; with this form just make up a list just consisting of the form.
		(or (assoc form1 (cdr p))
		    (list form1))
		(list form1 ))
;; Otherwise if a form was not specified, look up how many entries
;; in the operator dictionary for this operator.
	    (let* ((dict-entries (if p
				     (length (cdr  p))
				     0)))
	      (if (= 0 dict-entries)
;; If there are none, the operator is infix.
		  (list "infix" )
		  (if (= 1 dict-entries)
;; If there is one, return that.
		      (car (cdr p))
;; If there is more than one, look where we are in the parent list
;; to decide which one to use.
		      (let ((d (assoc (if (string=? "MROW" (gi pnt))
					  (let ((l (node-list-length (children pnt)))
						(n (absolute-child-number (current-node))))
					    (if (> l 1)
						(if (= n 1)
						    "prefix"
						    (if (= n l)
							"postfix"
							"infix"))
						"infix"))
					  "infix")
				      (cdr p))))
			(if d
			    d
			    (car (cdr p))))
		      )))	; let*
	    ))		; form
;; Fence, from an attribute or out of the dictionary
       (fence1 (attribute-string "fence"))
       (fence (math-my-debug (if fence1
			    (string=? fence1 "fence-true")
			    (cdr (or (assoc 'fence (cdr form))
				     '( #t . #f))))))
       )				; end of let* settings
;; Having done all that work, ignore all the attributes and just process
;; the children.
    (make  math-sequence
      (process-children-trim)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Helper function that just sticks a number of copies of a sosofo
;;; on to the current sosofo. Used to pad out grids.

(define (repeat-sosofo n sos)
  (if (equal? n 0)
      (empty-sosofo)
      (sosofo-append sos (repeat-sosofo (- n 1) sos))))

;;; mtable
;;; Spanning cells not supported as I need to use dsssl grid flow objects
;;; which don't span. Dsssl table flow objects do support spanning but they
;;; are display objects and can't be inlined so are no good for this.

(element mml:mtable
;; Preliminary pass through the table to count the number of columns
  (let ((cols (node-list-reduce
			 (children(current-node))
			 (lambda (a b)
			   (if (string=?  "MTR" (gi b))
			       (max a (node-list-length (children b)))
			       a))
			 0)))
    (make grid
      grid-n-columns: cols
;; Now main pass, making grid-cell flow object specifications.
      (node-list-reduce
       (children(current-node))
       (lambda (a b)
	 (if (string=?  "MTR" (gi b))
             (sosofo-append
	      a
	      (node-list-reduce 
	       (children b)
	       (lambda (x y)
		 (sosofo-append
		  x
		  (make grid-cell (process-node-list y))))
	       (empty-sosofo))
	      (repeat-sosofo
	       (- cols (node-list-length (children b)))
	       (make grid-cell (empty-sosofo))))
             (sosofo-append
	      a
	      (make grid-cell (process-node-list b))
	      (repeat-sosofo (- cols 1) (make grid-cell (empty-sosofo))))))
       (empty-sosofo)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; semantics

;;; Helper function to dig out a MathML-Presentation child
;;; if it exists.

(define (get-presentation nl)
  (if (node-list-empty? nl)
    #f
  (if (and 
       (string=?  "ANNOTATION-XML" (gi (node-list-first nl)))
       (string=?  "MATHML-PRESENTATION" 
		  (or (attribute-string "encoding" (node-list-first nl)) "")))
      (children (node-list-first nl))
      (get-presentation (node-list-rest nl)))))

;;; Typeset the body and ignore all annotations, unless there
;;; is an annotation-xml giving  MathML-Presentation 
;;; in which case use that and ignore everything else.

(element semantics
  (let* ((nl  (children (current-node))))
    (process-node-list (or 
			(get-presentation (node-list-rest nl))
			(node-list-first nl)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Content MathML

;;; In contrast to Presentation forms, Content MathML has vary few
;;; elements defined in the dsssl, instead the apply element explitly
;;; processes its children and calls [allegedly] suitable functions.

;;; It would be nice to use the dsssl transformation language to explicitly
;;; map content forms to presentation,
;;; but the style language isn't available, so... 


;;; Helper function dealing with bound variables on sums and products
;;; Checks for lowlimit child.

(define (dobvar product-char-sosofo equal-char-sosofo b r)
  (let ((first (gi (node-list-first r))))
    (cond ((string=? "lowlimit" first)
	   (bvarl  product-char-sosofo
	    (sosofo-append
	     b
	     equal-char-sosofo
	     (process-node-list (node-list-first r)))
	   (node-list-rest r)))
	   ((string=? "condition" first)
           (bvarl
	    product-char-sosofo
	    (process-node-list (node-list-first r))
	    (node-list-rest r)))
	   (else (bvarl product-char-sosofo (empty-sosofo) r)))))


;; Second function, checks for uplimit

(define (bvarl product-char-sosofo low r)
    (if (string=? "uplimit"  (gi (node-list-first r)))
	(bvaru 
	 product-char-sosofo
	 low
	 (process-node-list (node-list-first r))
	 (node-list-rest r))
	(bvaru
	 product-char-sosofo
	 low
	 (empty-sosofo)
	 r)))

;; third function, add the found limits to the sum or product character

<!--***changed by Herta: old function ***-->

<!--(define (bvaru product-char-sosofo low high r)
 (sosofo-append
 (make script
  (make math-sequence label: 'post-sub low)
  (make math-sequence label: 'post-sup high)
  product-char-sosofo)
    (process-node-list r))-->


<!--***changed by herta: new function ***-->
(define (bvaru product-char-sosofo low high r)
 (sosofo-append
 (make script
  (make math-sequence label: 'mid-sub low)
  (make math-sequence label: 'mid-sup high)
  (literal product-char-sosofo))
 
 (if (equal? "Max" product-char-sosofo)
    (make fence
     (make math-sequence label: 'open
	   (literal "("))
     (make math-sequence label: 'close
	   (literal ")"))
     (make math-sequence
       (process-node-list r)))
     (process-node-list r))))



;;; General function to typeset a function with braced arguments.

(define (apply-braced-fn fn args)
 (make math-sequence
  (sosofo-append
   (literal fn)
   (make fence
     (make math-sequence label: 'open
	   (literal "("))
     (make math-sequence label: 'close
	   (literal ")"))
     (make math-sequence 
       (process-node-list args))))))


;;; product/ sum/
;;; Look for bound vars, then call above helper function.

(define (apply-product product-char-sosofo equal-char-sosofo args)
  (make math-sequence
    (let ((b (node-list-first args)))
      (if (string=? (gi b) "bvar")
          (dobvar
	   product-char-sosofo
	   equal-char-sosofo
	   (process-node-list b)
	   (node-list-rest args))
	  (sosofo-append
	   product-char-sosofo
	   (process-node-list b))))))

;;; log/
;;; Just needs to check for logbase/

(define (apply-log args)
  (make math-sequence
    (if (equal? "LOGBASE" (gi (node-list-first args)))
	(sosofo-append
	 (make script
	   (sosofo-append
	    (literal "log")
	    (make math-sequence
	      label: 'post-sub
	      (process-node-list (node-list-first args)))))
	 (apply-braced-fn "" (node-list-rest args)))
	(apply-braced-fn "log" args))))

;;; root/ sqrt/
;;; Checks for degree/

(define (apply-root args)
  (make radical
    (if (equal? "DEGREE" (gi (node-list-first args)))
     (sosofo-append
      (make math-sequence
         (process-node-list (node-list-rest args)))
      (make math-sequence
        label: 'degree
         (process-node-list (node-list-first args))))
     (make math-sequence
       (process-node-list args)))))


;;; forall/ exists/
;;; Typeset any bound variables and conditions,
;;; separated by commas, then do a full stop and finally
;;; the body.

(define (apply-forall op c)
(make math-sequence
 (sosofo-append
   op
  (process-node-list (node-list-first c))
  (node-list-reduce
   (node-list-rest c)
   (lambda (result n)
     (sosofo-append
      result
      (literal
       (if (or (equal? "bvar" (gi n)) (equal? "condition" (gi n)))
       ":"																	<!--***changed by Herta***-->
       "|"))																<!--***changed by Herta***-->
      (process-node-list n)))
   (empty-sosofo)))))
     

;;; int/

<!--***changed by Herta: old function ***-->
<!--(define (apply-int  args)
  (make math-sequence
    (let ((b (node-list-first args)))
      (if (string=? (gi b) "bvar")
          (dobvar
	      (literal "\integral")
	   (empty-sosofo)
	   (empty-sosofo)
	   (node-list-rest args))
	  (sosofo-append
	   (literal "\integral")
	   (process-node-list b))))))-->

<!-- *** new function ***-->
(define (apply-int  args)
  (make math-sequence
    (let ((b (node-list-first args)))
      (if (string=? (gi b) "bvar")
      (make sequence
         (dobvar
	        "\integral"
	        (empty-sosofo)
	        (empty-sosofo)
	        (node-list-rest args))
         (literal "d")
         (process-node-list b))
	  (sosofo-append
	   (literal "\integral")
	   (process-node-list b))))))


;;; power/                     <!--***changed by Herta***-->

(define (apply-power args)
  (make script
     (sosofo-append
      (make math-sequence
      (if (and (equal? 1 (node-list-length (node-list-first args)))
		    (or (equal?  "ci" (gi (node-list-first args)))
              (equal?  "cn" (gi (node-list-first args)))))
            (process-node-list (node-list-first args))
      		(make fence
	           (make math-sequence label: 'open
		          (literal "("))
	           (make math-sequence label: 'close
		          (literal ")"))
	           (make math-sequence 
		 		    (process-node-list (node-list-first args)))))
      (make math-sequence
        label: 'post-sup
         (process-node-list (node-list-rest args)))))))


;;; factorial/                     <!--***added by Herta***-->

(define (apply-factorial args)
     (sosofo-append
      (make math-sequence
      (if (and (equal? 1 (node-list-length (node-list-first args)))
		    (or (equal?  "ci" (gi (node-list-first args)))
              (equal?  "cn" (gi (node-list-first args)))))
            (process-node-list (node-list-first args))
      		(make fence
	           (make math-sequence label: 'open
		          (literal "("))
	           (make math-sequence label: 'close
		          (literal ")"))
	           (make math-sequence 
		 		    (process-node-list (node-list-first args)))))
      (literal "!"))))



;;; inverse/

(define (apply-inverse fname args)
  (make script
     (sosofo-append
      (make math-sequence
         (if (and (equal? 1 (node-list-length (node-list-first args)))
		    (or (equal?  "ci" (gi (node-list-first args)))
              (equal?  "cn" (gi (node-list-first args)))))
            (process-node-list (node-list-first args))
      		(make fence
	           (make math-sequence label: 'open
		          (literal "("))
	           (make math-sequence label: 'close
		          (literal ")"))
	           (make math-sequence 
		 		    (process-node-list (node-list-first args)))))
      (make math-sequence
        label: 'post-sup
	     (case fname
	        (("inverse") (literal "-1"))
	        (("transpose") (literal "t"))))))))

;;; over bars (don't work)
(define (apply-over args)
  (make math-sequence
       (process-node-list args))) ;; fix!

;;; max/ and min/

(define (apply-max product-char-sosofo equal-char-sosofo args)
 (make math-sequence
    (let ((b (node-list-first args)))
      (if (string=? (gi b) "bvar")
          (dobvar
	   product-char-sosofo
	   equal-char-sosofo
	   (process-node-list b)
	   (node-list-rest args))
	  (sosofo-append
	   (literal product-char-sosofo)
	   (make fence
		  (make math-sequence label: 'open
	       (literal "("))             
        (make math-sequence label: 'close
	       (literal ")"))
        (make math-sequence
          (process-node-list b))))))))




<!--  (make math-sequence
    (sosofo-append
     (make unmath (literal fname))
     (do-set "(" ")" args))))-->                 ;;***changed by Herta new parentheses!! ***

;;; gcd/

(define (apply-gcd fname args)
  (make math-sequence
    (sosofo-append
     (make unmath (literal fname))
     (do-set "(" ")" args))))

;;; abs/

(define (apply-abs args)
 (make fence
   (make math-sequence label: 'open
	 (literal "|"))
   (make math-sequence label: 'close
	 (literal "|"))
   (make math-sequence
     (process-node-list args))))


;;; norm/        --*** added by Herta ***--

(define (apply-norm args)
 (if (equal? "DEGREE" (gi (node-list-first args)))
  (make script
   (sosofo-append
    (make math-sequence
     (make fence
       (make math-sequence label: 'open
 	     (literal "\U-2016"))             ;;***muß normalerweise "\U-01C1" heißen - geht aber noch nich
       (make math-sequence label: 'close 
     	  (literal "\U-2016"))
       (make math-sequence
         (process-node-list (node-list-rest args)))))
     (make math-sequence
        label: 'post-sub
	     (process-node-list (node-list-first args)))))

  (make fence
   (make math-sequence label: 'open
	  (literal "\U-2016"))             ;;***muß normalerweise "\U-01C1" heißen - geht aber noch nich
   (make math-sequence label: 'close
	 (literal "\U-2016"))
   (make math-sequence
     (process-node-list args)))))


;;; diff/
;;; Trick is to get hold of the degree. 
;;; Recommendation is for 
;;;  d f
;;;  --- (x)
;;;  d x
;;;
;;; but I don't currently test for special functions which are
;;; just a single identifier like that, so currently you get
;;;  d              d^3       
;;;  --- f(x)  and  ---   f(x)
;;;  d x     	    d x^3     

(define (apply-diff args)
  (sosofo-append
  (let*  ((f (node-list-first args)))
    (if (equal? "bvar" (gi f))
	(let ((d (node-list-filter
		  (lambda (n) (equal? "DEGREE" (gi n)))
                  (children f))))
	  (make fraction
	    (sosofo-append
	     (make math-sequence
	       label: 'numerator
	       (make script 
		 (sosofo-append
		  (literal "d")
		  (make math-sequence
		    label: 'post-sup 
		    (process-node-list d))))))
	    (make math-sequence
	      label: 'denominator
	      (sosofo-append
	       (literal "d")
               (make script
		 (sosofo-append
		  (process-node-list
		   (node-list-filter
		    (lambda (n) (not (equal? "DEGREE" (gi n))))
		    (children f)))
		  (make math-sequence
		    label: 'post-sup 
		    (process-node-list d))))))))
	(make fraction
	  (sosofo-append
	   (make math-sequence
	     label: 'numerator
	     (literal "d"))
	   (make math-sequence
	     label: 'denominator
	     (sosofo-append
	      (literal "d")
	      (process-node-list f)))))))
  (process-node-list (node-list-rest args))))


;;; partial differentiation
;;;
;;; Trick here is to amalgamate the degrees
;;; in the numerator:
;;;       d^3       
;;;  ---------- f(x)
;;;  d^2 x  d y  
;;;

;; helper function that sets the bit of the denominator
;;  corresponding to one bound variable.

(define (do-one-partial nl)
  (let ((d (node-list-filter
			(lambda (n) (equal? "DEGREE" (gi n)))
			(children nl))))
    (make math-sequence
      label: 'denominator
      (sosofo-append
       (literal "d")
       (make script
	 (sosofo-append
	  (process-node-list
	   (node-list-filter
	    (lambda (n) (not (equal? "DEGREE" (gi n))))
	    (children nl)))
	  (make math-sequence
	    label: 'post-sup 
	    (process-node-list d))))))))


;;next function added by Herta for patial differentiation


(define (do-one-partial-diff nl)
  (let ((d (node-list-filter
			(lambda (n) (equal? "DEGREE" (gi n)))
			(children nl))))
    (make math-sequence
      label: 'denominator
      (sosofo-append
       (literal "\partial-differential")
       (make script
	 (sosofo-append
	  (process-node-list
	   (node-list-filter
	    (lambda (n) (not (equal? "DEGREE" (gi n))))
	    (children nl)))
	  (make math-sequence
	    label: 'post-sup 
	    (process-node-list d))))))))


;;; partialdiff
;;; I suspect that I could make use of node-list-filter-by-gi
;;; from dblib here, but I only just noticed that function...

(define (apply-partialdiff args)
  (sosofo-append
   (let*  ((bvars
	    (node-list-filter
	     (lambda (n) (equal? "bvar" (gi n)))
	     args))
;; totald is the sum of all the degrees [with 1 being supplied
;; as a default in each case] There may be an easier way to get that.
	  (totald
	   (node-list-reduce
	    bvars
	    (lambda (result n)
	      (+ result
		 (or
		  (string->number 
		   (string-nospace
		    (data 
		     (node-list-filter
		      (lambda (nn) (equal? "DEGREE" (gi nn)))
		      (children n)))))
                     1)))
	    0)))
;; Now it is easy, set a fraction with the partial and the bound vars
;; then put out the body.
     (make fraction
       (sosofo-append
	(make math-sequence
	  label: 'numerator
	  (make script 
	    (sosofo-append
	     (literal "\partial-differential")
		  (if (equal? totald 1)
			  (empty-sosofo)
   	     (make math-sequence
	          label: 'post-sup 
	          (literal (number->string totald)))))))
	(make math-sequence
	  label: 'denominator
	  (node-list-reduce
	   bvars
	   (lambda (result n)
	     (sosofo-append
;; something is reversing the args, so I'll reverse them back
	      (do-one-partial-diff n)
	      result))
	   (empty-sosofo))))))
   (process-node-list (node-list-filter
		       (lambda (n) (not (equal? "bvar" (gi n))))
		       args))))



;;exp    															<!--***added by Herta***-->
(define (apply-exp args)
  (make script
     (sosofo-append
      (make math-sequence
      	(literal "e")
	      (make math-sequence
        		label: 'post-sup
         	(process-node-list args))))))


;;determinant
(define (apply-determinant args)
   (make math-sequence
     (literal "det")
     (make fence
	       (make math-sequence label: 'open
		     (literal "("))
	       (make math-sequence label: 'close
		     (literal ")"))
	       (make math-sequence 
		      (process-node-list args)))))




;;; apply
;;; This is the main switch for Content MathML
;;; Mainly just a case statement on th egi of the first child
;;; with some fall back code in case the function argument isn't known.
(element apply 
 (let* ((nl (children (current-node)))
	(f (node-list-first nl))
        (fname (gi f))
        (args (node-list-rest nl)))
   (case fname
     (("product")
      (apply-product "\n-ary-product" (literal "=") args))          ;;changed by Herta: old argument: (literal "\n-ary-product")
     (("sum")
      (apply-product "\n-ary-summation" (literal "=") args))         ;;changed by Herta: old argument: (literal "\n-ary-summation")
     (("limit")
      (apply-product "lim" (literal "\rightwards-arrow") args))       ;;changed by Herta: old argument: (literal "lim")
     (("int")
      (apply-int  args))
     (("inverse" "transpose")
      (apply-inverse fname  args))
     (("power")
       (apply-power  args))                  <!-- if-->
     (("exp")
       (apply-exp args))
     (("forall")
      (apply-forall (literal "\for-all")  args))
     (("exists")
      (apply-forall (literal "\there-exists")  args))
     (("plus" "minus" "union" "intersect" "compose" "divide" "rem" "and" "or")          ;;***"'AND' 'OR'"  added by Herta***
      (do-nary-binop fname  args))
     (("times")
      (if (equal? "1" (attribute-string "visual"))
          (do-nary-binop "timesvis" args)
          (do-nary-binop fname args)))
	  (("fn")
		(make sequence 
		   (if (equal? "0" (attribute-string "visual"))
		    (process-children)
          (sosofo-append
				(process-node-list f)
				(make fence
	          (make math-sequence label: 'open
		        (literal "("))
	          (make math-sequence label: 'close
		        (literal ")"))
	          (make math-sequence 
		        (process-node-list args)))))))                                   ;;HACK by Herta	
     (("mean" "conjugate")
      (apply-over args))
     (("factorial")
      (apply-factorial args))                                ;;changed by Herta
     (("root")
      (apply-root args))
     (("sdev")
      (apply-braced-fn "\greek-small-letter-sigma" args))
     (("median" "mode")
      (apply-braced-fn fname args))
     (("max")                                                ;;changed by Herta: old argument: (literal "Max")
      (apply-max "Max" (literal "=") args))
     (("min")
      (apply-max "Min" (literal "=") args))         ;;changed by Herta: old argument: (literal "Min")
     (("gcd")
      (apply-gcd fname args))
     (("log")
      (apply-log args))
     (("abs")
      (apply-abs  args))
     (("norm")                                   ;;***added by Herta***
      (apply-norm  args))
     (("diff")
      (apply-diff  args))
     (("partialdiff")
      (apply-partialdiff args))
     (("determinant")                             ;;***added by Herta***
      (apply-determinant args))          
     (("var")
      (make script
      (make math-sequence
        label: 'post-sup (literal "2"))
        (apply-braced-fn "\greek-small-letter-sigma" args)))
     (("transpose")
      (make script
      (make math-sequence
        label: 'post-sup (literal "t"))
        (process-node-list args)))
     (else 
      (make math-sequence
	(sosofo-append
	 (process-node-list f)
	 (if (and (equal? 
		   1 (node-list-length args))
		  (equal?  "ci" (gi args)))
	     (make math-sequence (process-node-list args))
	     (make fence
	       (make math-sequence label: 'open
		     (literal "("))
	       (make math-sequence label: 'close
		     (literal ")"))
	       (make math-sequence 
		 (process-node-list args))))))))))


;;; reln
;;; relations are similar to apply.
;;; Currently there is a double test as this top level
;;; function bunches all binops together, then  do-binary-reln          
;;; tests again to separate them, perhaps not the best idea.

(element reln
  (let* ((r (gi (node-list-first (children (current-node)))))
	 (c (node-list-rest(children (current-node)))))
  (case r
   (("neq" "implies" "in" "notin" "notsubset" "notprsubset" "tendsto")
     (do-binary-reln r c))
    (("eq" "leq" "lt" "geq" "gt" "subset" "prsubset")
    (do-nary-reln r c))
   )))


;;; Typeset a binary infix relation
  
(define (do-binary-reln r c)
 (sosofo-append
  (process-node-list (node-list-first  c))
  (literal (case r
	     (("neq") "\not-equal-to")
	     (("implies") "\rightwards-double-arrow")          <!--***changed by Herta***-->
	     (("in") "\element-of")
	     (("notin") "\not-an-element-of")
	     (("notsubset") "\not-a-subset-of")
	     (("notprsubset") "PR\not-a-subset-of"); fix!
	     (("tendsto") "\rightwards-arrow"))); fix for type attribute
  (process-node-list (node-list-rest  c))))

;;; Typeset a nary infix relation
;;; repeating the operator as many times as needed.

(define (do-nary-reln r c)
 (let ((rs (literal 
	    (case r
	      (("eq") "\equals-sign")
	      (("leq") "\less-than-or-equal-to")
	      (("lt") "\less-than-sign")
	      (("geq") "\greater-than-or-equal-to")
	      (("gt") "\greater-than-sign")
	      (("subset") "\subset-of")
	      (("prsubset") "PR\subset-of")))));fix!
   (sosofo-append
    (process-node-list (node-list-first c))
    (node-list-reduce
     (node-list-rest c)
     (lambda (result n)
       (sosofo-append
	result
	rs
	(process-node-list n)))
     (empty-sosofo)))))
     
;;; typeset an nary operator, as for relations.

(define (do-nary-binop op c)
 
 (let ((ops (literal 
	    (case op
	      ((",") ",")
	      (("plus") "\plus-sign")
	      (("timesvis") "\U-22C5")
         (("times") "")                                    <!--***changed by Herta - has to be fixed again!!***-->
	      (("union") "\union")
	      (("intersect") "\intersection")
	      (("minus") "\minus-sign")                         
	      (("divide") "/")                    <!--***changed by Herta - was : (("divide") "\division-slash")***-->
	      (("rem") "mod") ;fix                        
         (("and") "\U-2227")                                      <!--***added by Herta***-->
         (("or") "\U-2228")                                       <!--***added by Herta***-->
	      (("compose") "\U-2218") ))))
 (sosofo-append
  (let* ((nl (children (node-list-first c)))                                       <!--***all changed by Herta***-->
	(f (node-list-first nl))
        (fname (gi f))
        (args (node-list-rest nl)))
        (if (and (or (equal? "and" fname) (equal? "or" fname) (equal? "times" op) (equal? "timesvis" op) (equal? "divide" op) (equal? "rem" op) (equal? "compose" op))   <!--***UNION, intERSECT?***-->
                 (or (equal? "PLUS" fname) (equal? "minus" fname) (equal? "or" fname) (equal? "and" fname)))
              (make fence
	           	(make math-sequence label: 'open
		           (literal "("))
	            (make math-sequence label: 'close
		           (literal ")"))
	            (make math-sequence
		 		     (process-node-list (node-list-first c))))
              (process-node-list (node-list-first c))))
  (node-list-reduce
   (node-list-rest c)
   (lambda (result n)
     (sosofo-append
      result
      ops
      (let* ((nl (children (node-list-first n)))
	     (f (node-list-first nl))
          (fname (gi f))
          (args (node-list-rest nl)))
		    (if (or (and (or (equal? "and" fname) (equal? "or" fname) (equal? "times" op) 
                           (equal? "timesvis" op) (equal? "divide" op) (equal? "rem" op) (equal? "compose" op))  
                       (or (equal? "and" fname) (equal? "or" fname) (equal? "PLUS" fname) (equal? "minus" fname)))                               <!--***UNION, intERSECT?***-->
                  (and (equal? "minus" op)
                       (or (equal? "plus" fname) (equal? "minus" fname))))
                 (make fence
	           	     (make math-sequence label: 'open
		                (literal "("))
	                 (make math-sequence label: 'close
		                (literal ")"))
	                 (make math-sequence 
		 		          (process-node-list n)))
                  (process-node-list n)))))                                      <!--***changed by Herta end***-->
   (empty-sosofo)))))
     



;;;;;;;;;;;;;;;;

;;;; ident

(element ident
  (make math-sequence
    (literal "id")))



;;; trig
;;; Currently these done with th eelement declarations
;;; and the fallback case of apply typesetting the arguments.
;;; may need to change that.

(define (do-sin sin) 
  (make unmath 
    font-posture: 'upright
    (literal sin)))

(element sin
  (do-sin "sin"))

(element cos
  (do-sin "cos"))

(element tan
  (do-sin "tan"))

(element sec
  (do-sin "sec"))

(element csc
  (do-sin "csc"))

(element cot
  (do-sin "cot"))

(element sinh
  (do-sin "sinh"))

(element cosh
  (do-sin "cosh"))

(element tanh
  (do-sin "tanh"))

(element sech
  (do-sin "sech"))

(element csch
  (do-sin "csch"))

(element coth
  (do-sin "coth"))

(element arcsin
  (do-sin "arcsin"))

(element arccos
  (do-sin "arccos"))

(element arctan
  (do-sin "arctan"))

;;;;

;;; natural log and determinant are honourary trig functions.
;;; I should make det applied to an explicit table use | | rather
;;; than det ( ) I think.

(element ln
  (do-sin "ln"))


;;(element determinant
;;  (do-sin "det"))



;;;;;;;;;;;

;;; sets and lists

(element mml:mset
 (do-set "{" "}" (children (current-node))))

(element list
 (do-set "[" "]" (children (current-node))))

(element vector
 (do-set "(" ")" (children (current-node))))

(element interval
  (case (attribute-string "closure" (current-node))
    (("closed")
     (do-set "[" "]" (children (current-node))))
    (("open")
     (do-set "(" ")" (children (current-node))))
    (("open-closed")
     (do-set "(" "]" (children (current-node))))
    (("closed-open")
     (do-set "[" ")" (children (current-node))))))


;;; helper function for sets.
;;; Checks for two styles, comma separated explicit,
;;; or via conditions and bound vars.

(define (do-set left right args)
 (make fence
   (make math-sequence label: 'open
	 (literal left))
   (make math-sequence label: 'close
	 (literal right))
   (make math-sequence
     (if (equal? "bvar" (gi (node-list-first args)))
	 (if (equal? 2 (node-list-length  args))
;; typeset bvar
	     (sosofo-append
	      (process-node-list (children (node-list-first args)))
	      (literal "|")
	      (process-node-list (children (node-list-rest args))))
;; else dont
	     (sosofo-append
	      (process-node-list (node-list-rest(node-list-rest args)))
	      (literal "|")
	      (process-node-list (children (node-list-first (node-list-rest args))))))
;; else stick in some commas
	 (do-nary-binop "," args)	 ))))


;;;;;;;

;;; lambda [lamda according to Jade ?]

(element lambda
 (sosofo-append
  (make math-sequence (literal "\greek-small-letter-lamda"))
   (do-nary-binop ","
    (node-list-filter
     (lambda (n) (equal? "bvar" (gi n)))
     (children (current-node))))
   (literal ".")
   (process-node-list
    (node-list-filter
     (lambda (n) (not (equal? "bvar" (gi n))))
     (children (current-node))))))


;;;;;;;;;;;;;;;

;;; cn

;;; what to do if you see a sep element for rationals..

(define (sep-rational top bottom)
  (make fraction
    (sosofo-append
     (make math-sequence
       label: 'numerator
       (process-node-list top))
     (make math-sequence
       label: 'denominator
       (process-node-list bottom)))))


;;; and for complex cartesian 

(define (sep-complex-cartesian top bottom)
  (make math-sequence
    (sosofo-append
     (make math-sequence
       (process-node-list top))
     (make math-sequence (literal "+"))
     (make math-sequence
       (process-node-list bottom))
     (make math-sequence (literal "i")))))


;;; and polar

(define (sep-complex-polar top bottom)
  (make math-sequence
    (sosofo-append
     (make math-sequence (literal "Polar"))
     (make fence
       (sosofo-append
      (make math-sequence
         label: 'open
          (literal "("))
      (make math-sequence
         label: 'close
          (literal ")"))
       (sosofo-append
	(make math-sequence
	  (process-node-list top))
	(make math-sequence (literal ","))
	(make math-sequence
	  (process-node-list bottom))))))))

;;; Function to collect node list before and after sep element
;;; then finally call one of the above functions to typeset the
;;; two halves of the number.

(define (do-sep fn flag top bottom args)
  (if (equal? 0 (node-list-length args)); why can't I use null
      (fn top bottom)
      (if (equal? "sep" (gi (node-list-first args)))
	  (do-sep fn #f top bottom (node-list-rest args))
          (if flag
	      (do-sep fn flag
		      (node-list
		       top
		       (node-list-first args))
		      bottom
		      (node-list-rest args))
	      (do-sep fn flag
		      top
		      (node-list
		       bottom
		       (node-list-first args))
		      (node-list-rest args))))))


;;; cn
;;; subscript with the base unless it is 10, or call one of the above functions
;;; to start looking for sep element.

(element cn
  (case (or (attribute-string "type" (current-node))
            "real")
    (("real" "constant" "vector" "matrix")
     (make math-sequence (process-node-list (children (current-node)))))

    (("integer")
     (if (equal? "10" (attribute-string "base" (current-node)))
	 (make math-sequence
	   (process-node-list (children (current-node))))	
	 (make script
	   (sosofo-append
	    (make math-sequence
	      (process-node-list (children (current-node))))
	    (make math-sequence
	      label: 'post-sub
	      (literal 
	       (attribute-string "base" (current-node))))))))
    (("rational")
     (do-sep sep-rational
	     #t
	     (empty-node-list)
	     (empty-node-list)
	     (children (current-node))))
    (("complex-cartesian")
     (do-sep sep-complex-cartesian
	     #t
	     (empty-node-list)
	     (empty-node-list)
	     (children
	      (current-node))))
    (("complex-polar")
     (do-sep sep-complex-polar
	     #t
	     (empty-node-list)
	     (empty-node-list)
	     (children (current-node))))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; ought to share code with mtable, probably

(element mml:matrix
  (let* ((c (math-my-debug (children (current-node))))
	 (n (math-my-debug  (node-list-length (children (node-list-first c))))))
    (make fence
      (make math-sequence label: 'open
	    (literal "("))
      (make math-sequence label: 'close
	    (literal ")"))
      (make grid
	grid-n-columns: n
	(node-list-reduce
	 c
	 (lambda (a b)
	   (sosofo-append
	    a
	    (node-list-reduce
	     (children b)
	     (lambda (x y)
	       (sosofo-append
		x
		(make grid-cell
		  (process-node-list y))))
	     (empty-sosofo))))
	 (empty-sosofo))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Always displayed, so far.
;;; Should do inline as well.

;;for equation only

(element (equation mml:math)
 (make display-group
   min-leading: 2pt
 (make math-sequence
   math-display-mode: 'display
     (process-children-trim))))



(element (informalequation mml:math)
 (make display-group
   min-leading: 2pt
 (make math-sequence
   math-display-mode: 'display
     (process-children-trim))))

(element (inlineequation mml:math)
 (make sequence
 (make math-sequence
   math-display-mode: 'inline
     (process-children-trim))))

;;ohne gewaehr, da cut and paste! attribute ueberpruefen !!!
;; der typ hat ein haufen elemente vergessen
(element fn (process-children-trim))
(element uplimit  (process-children-trim))
(element lowlimit (process-children-trim))
(element bvar  (process-children-trim))
(element degree (process-children-trim))
;(element forall (process-children-trim))
;(element exists  (process-children-trim))



<!-- ***following function changed by Herta *** -->

(element condition  
  (if (equal? 1 (node-list-length (children (current-node))))
      (process-children-trim)
      (make math-sequence
         (process-node-list (node-list-first (children (current-node))))
         (literal " : ")
         (process-node-list (node-list-rest (children (current-node)))))))


(element logbase (process-children-trim))

(element ci
  (make math-sequence
    font-posture: 
    (let ((fnt 
	   (if (attribute-string "fontstyle")
	      (attribute-string "fontstyle")
	      "italic")))
      (if (equal? "style-normal" fnt)
	  'upright
	  (if (equal? "italic" fnt)
	      'italic 
	      #f)))
    (process-children-trim)

))



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define math-my-debug (lambda (x) x))

<!--
(define math-my-debug2 
 (lambda (x)
 ((lambda (a b) b)
 (debug (debug-body x))
 x)))


(define max-node-list-debug-length 5)

(define (debug-body x)
  `(,(cond ((node-list? x)
	    (if (node-list-empty? x)
		'empty-node-list
		`( ,(if (named-node-list? x)
			'named-node-list
			'node-list)
		   ,(node-list-length x)
		   ,(node-list-reduce 
		      (node-list-head x max-node-list-debug-length) 
		      (lambda (result n) 
			(string-append result
				       (cond ((gi n) (string-append "<" (gi n) ">" ))
					     ((equal? 'data-char (node-property 'class-name n))  (data n))
					     (else "<?>"))))
		      "" ))))
	   ((sosofo? x) 'sosofo)
	   ((procedure? x) 'procedure)
	   ((style? x) 'style)
	   ((address? x) 'address)
	   ((color? x) 'color)
	   ((color-space? x) 'color-space)
	   ((display-space? x) 'display-space)
	   ((inline-space? x) 'inline-space)
	   ((glyph-id? x) 'glyph-id)
	   ((glyph-subst-table? x) 'glyph-subst-table)
	   ((boolean? x) 'boolean)
	   ((symbol? x) 'symbol)
	   ((list? x) 'list)
	   ((pair? x) 'pair)
	   ((char? x) 'char)
	   ((string? x) 'string)
	   ((quantity? x) 'quantity)
	   ((keyword? x) 'keyword)
	   (else 'other))
    ,x))

-->
<!-- Some bits from the dsssl report, mainly from the
     Mulberry site
  -->

(define (node-list-filter proc nl)
  (node-list-reduce nl
                    (lambda (result snl)
                      (if (proc snl)
                          (node-list snl result)
                          result))
                    (empty-node-list)))



(define (node-list-head nl i)
  (if (zero? i) 
      (empty-node-list)
      (node-list (node-list-first nl)
;;; page 136 of dsssl spec appears to be wrong...
		 (node-list-head (node-list-rest nl)
				 (- i 1)))))

(define (zero? x) (equal? x 0))

(define (attribute name nl) 
 (node-list-map (lambda (snl) (named-node name (attributes snl))) nl))

(define (attribute-value name nl)
      (node-list-property 'value (attribute name nl)))

(define (node-list-property prop nl) (node-list-map (lambda (snl)
 (node-property prop snl default: (empty-node-list))) nl))

