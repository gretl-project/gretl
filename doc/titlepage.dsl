;; BOOK titlepage

(mode book-titlepage-recto-mode
  (element abbrev
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element abstract
    (make display-group
      use: book-titlepage-recto-style
      quadding: 'start
      ($semiformal-object$)))

  (element (abstract title) (empty-sosofo))

  (element (abstract para)
    (make paragraph
      use: book-titlepage-recto-style
      quadding: 'start
      (process-children)))

  (element address 
    (make display-group
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (with-mode titlepage-address-mode 
	($linespecific-display$ %indent-address-lines% %number-address-lines%))))

  (element affiliation
    (make display-group
      use: book-titlepage-recto-style
      (process-children)))

  (element author
    (let ((author-name  (author-string))
	  (author-affil (select-elements (children (current-node)) 
					 (normalize "affiliation"))))
      (make sequence      
	(make paragraph
	  use: book-titlepage-recto-style
	  font-size: (HSIZE 3)
	  line-spacing: (* (HSIZE 3) %line-spacing-factor%)
	  space-before: (* (HSIZE 2) %head-before-factor%)
	  quadding: %division-title-quadding%
	  keep-with-next?: #t
	  (literal author-name))
	(process-node-list author-affil))))

  (element authorblurb
    (make display-group
      use: book-titlepage-recto-style
      quadding: 'start
      (process-children)))

  (element (authorblurb para)
    (make paragraph
      use: book-titlepage-recto-style
      quadding: 'start
      (process-children)))

  (element authorgroup
    (make display-group
      (process-children)))

  (element authorinitials
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element copyright
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (literal (gentext-element-name (current-node)))
      (literal "\no-break-space;")
      (literal (dingbat "copyright"))
      (literal "\no-break-space;")
      (process-children)))

  (element (copyright year)
    (make sequence
      (process-children)
      (if (not (last-sibling? (current-node)))
	  (literal ", ")
	  (literal (string-append " " (gentext-by) " ")))))
  
  (element (copyright holder) ($charseq$))

  (element date
    (make paragraph
      use: book-titlepage-recto-style
      space-before: (* (HSIZE 4) %head-before-factor%)
      quadding: %division-title-quadding%
      (process-children)))

  (element edition
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)
      (literal "\no-break-space;")
      (literal (gentext-element-name-space (gi (current-node))))))

  (element firstname
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element graphic
    (let* ((nd (current-node))
	   (fileref (attribute-string "fileref" nd))
	   (entityref (attribute-string "entityref" nd))
	   (format (attribute-string "format" nd))
	   (align (attribute-string "align" nd)))
      (if (or fileref entityref) 
	  (make external-graphic
	    notation-system-id: (if format format "")
	    entity-system-id: (if fileref 
				  (graphic-file fileref)
				  (if entityref 
				      (entity-generated-system-id entityref)
				      ""))
	    display?: #t
	    display-alignment: 'center)
	  (empty-sosofo))))

  (element honorific
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element keywordset
    (make paragraph
      quadding: 'start
      (make sequence
	font-weight: 'bold
	(literal "Keywords: "))
      (process-children)))

  (element (keyword)
    (make sequence
      (process-children)
      (if (not (last-sibling?))
	  (literal ", ")
	  (literal ""))))

  (element legalnotice
    (make display-group
      use: book-titlepage-recto-style
      ($semiformal-object$)))

  (element (legalnotice title) (empty-sosofo))

  (element (legalnotice para)
    (make paragraph
      use: book-titlepage-recto-style
      quadding: 'justify
      line-spacing: (* 0.8 (inherited-line-spacing))
      font-size: (* 0.8 (inherited-font-size))
      (process-children)))

  (element modespec (empty-sosofo))

  (element orgdiv
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element orgname
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element othercredit
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element othername
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element pagenums
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element printhistory
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element pubdate
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element publisher
    (make display-group
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element publishername
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element releaseinfo
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element revhistory
    (make sequence
      (make paragraph
	use: book-titlepage-recto-style
	space-before: (* (HSIZE 3) %head-before-factor%)
	space-after: (/ (* (HSIZE 1) %head-before-factor%) 2)
	(literal (gentext-element-name (current-node))))
      (make table
	before-row-border: #f
	(process-children))))
  
  (element (revhistory revision)
    (let ((revnumber (select-elements (descendants (current-node)) 
				      (normalize "revnumber")))
	  (revdate   (select-elements (descendants (current-node)) 
				      (normalize "date")))
	  (revauthor (select-elements (descendants (current-node))
				      (normalize "authorinitials")))
	  (revremark (select-elements (descendants (current-node))
				      (normalize "revremark"))))
      (make sequence
	(make table-row
	  (make table-cell
	    column-number: 1
	    n-columns-spanned: 1
	    n-rows-spanned: 1
	    start-indent: 0pt
	    (if (not (node-list-empty? revnumber))
		(make paragraph
		  use: book-titlepage-recto-style
		  font-size: %bf-size%
		  font-weight: 'medium
		  (literal (gentext-element-name-space (current-node)))
		  (process-node-list revnumber))
		(empty-sosofo)))
	  (make table-cell
	    column-number: 2
	    n-columns-spanned: 1
	    n-rows-spanned: 1
	    start-indent: 0pt
	    cell-before-column-margin: (if (equal? (print-backend) 'tex)
					   6pt
					   0pt)
	    (if (not (node-list-empty? revdate))
		(make paragraph
		  use: book-titlepage-recto-style
		  font-size: %bf-size%
		  font-weight: 'medium
		  (process-node-list revdate))
		(empty-sosofo)))
	  (make table-cell
	    column-number: 3
	    n-columns-spanned: 1
	    n-rows-spanned: 1
	    start-indent: 0pt
	    cell-before-column-margin: (if (equal? (print-backend) 'tex)
					   6pt
					   0pt)
	    (if (not (node-list-empty? revauthor))
		(make paragraph
		  use: book-titlepage-recto-style
		  font-size: %bf-size%
		  font-weight: 'medium
		  (literal (gentext-revised-by))
		  (process-node-list revauthor))
		(empty-sosofo))))
	(make table-row
	  cell-after-row-border: #f
	  (make table-cell
	    column-number: 1
	    n-columns-spanned: 3
	    n-rows-spanned: 1
	    start-indent: 0pt
	    (if (not (node-list-empty? revremark))
		(make paragraph
		  use: book-titlepage-recto-style
		  font-size: %bf-size%
		  font-weight: 'medium
		  space-after: (if (last-sibling?) 
				   0pt
				   (/ %block-sep% 2))
		  (process-node-list revremark))
		(empty-sosofo)))))))
  
  (element (revision revnumber) (process-children-trim))
  (element (revision date) (process-children-trim))
  (element (revision authorinitials) (process-children-trim))
  (element (revision revremark) (process-children-trim))

  (element shortaffil
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element subjectset (empty-sosofo))

  (element subtitle 
    (make paragraph
      use: book-titlepage-recto-style
      font-size: (HSIZE 4)
      line-spacing: (* (HSIZE 4) %line-spacing-factor%)
      space-before: (* (HSIZE 4) %head-before-factor%)
      quadding: %division-subtitle-quadding%
      keep-with-next?: #t
      (process-children-trim)))

  (element surname
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element title 
    (make paragraph
      use: book-titlepage-recto-style
      font-size: (HSIZE 8)
      line-spacing: (* (HSIZE 8) %line-spacing-factor%)
      space-before: (* (HSIZE 8) %head-before-factor%)
      space-after: (* (HSIZE 8) %line-spacing-factor%)
      quadding: %division-title-quadding%
      keep-with-next?: #t
      heading-level: (if %generate-heading-level% 1 0)
      (with-mode title-mode
	(process-children-trim))))

  (element titleabbrev (empty-sosofo))

  (element volumenum
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))
)

;; ARTICLE titlepage

(mode article-titlepage-recto-mode

  (element affiliation
    (make display-group
      use: book-titlepage-recto-style
      (process-children)))

  (element author
    (let ((author-name  (author-string))
	  (author-affil (select-elements (children (current-node)) 
					 (normalize "affiliation"))))
      (make sequence      
	(make paragraph
	  use: book-titlepage-recto-style
	  font-size: (HSIZE 3)
	  line-spacing: (* (HSIZE 3) %line-spacing-factor%)
	  space-before: (* (HSIZE 6) %head-before-factor%)
	  quadding: %division-title-quadding%
	  keep-with-next?: #t
	  (literal author-name))
	(process-node-list author-affil))))

  (element othercredit
    (let ((othercredit-name  (author-string))
	  (othercredit-contrib (select-elements (children (current-node)) 
					 (normalize "contrib"))))
      (make sequence      
	(make paragraph
	  use: book-titlepage-recto-style
	  font-size: (HSIZE 3)
	  line-spacing: (* (HSIZE 3) %line-spacing-factor%)
	  space-before: (* (HSIZE 6) %head-before-factor%)
	  quadding: %division-title-quadding%
	  keep-with-next?: #t
	  (literal othercredit-name))
	(process-node-list othercredit-contrib))))

  (element copyright
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (literal (gentext-element-name (current-node)))
      (literal "\no-break-space;")
      (literal (dingbat "copyright"))
      (literal "\no-break-space;")
      (process-children)))

  (element (copyright year)
    (make sequence
      (process-children)
      (if (not (last-sibling? (current-node)))
	  (literal ", ")
	  (literal (string-append " " (gentext-by) " ")))))
  
  (element (copyright holder) ($charseq$))

  (element date
    (make paragraph
      use: book-titlepage-recto-style
      space-before: (* (HSIZE 8) %head-before-factor%)
      quadding: %division-title-quadding%
      (process-children)))

  (element edition
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)
      (literal "\no-break-space;")
      (literal (gentext-element-name-space (gi (current-node))))))

;;  (element firstname
;;    (make paragraph
;;      use: book-titlepage-recto-style
;;      quadding: %division-title-quadding%
;;      (process-children)))

  (element legalnotice
    (make display-group
      use: book-titlepage-recto-style
      ($semiformal-object$)))

  (element (legalnotice title) (empty-sosofo))

  (element (legalnotice para)
    (make paragraph
      use: book-titlepage-recto-style
      quadding: 'justify
      line-spacing: (* 0.8 (inherited-line-spacing))
      font-size: (* 0.8 (inherited-font-size))
      space-before: (* 4 (inherited-line-spacing))
      (process-children)))

  (element modespec (empty-sosofo))

  (element orgdiv
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element orgname
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

;;  (element othercredit
;;    (make paragraph
;;      use: book-titlepage-recto-style
;;      font-size: (HSIZE 3)
;;      line-spacing: (* (HSIZE 3) %line-spacing-factor%)
;;      space-before: (* (HSIZE 6) %head-before-factor%)
;;      quadding: %division-title-quadding%
;;      (process-children)))

  (element othername
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element pubdate
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element releaseinfo
    (make paragraph
      use: book-titlepage-recto-style
      quadding: %division-title-quadding%
      (process-children)))

  (element subtitle 
    (make paragraph
      use: book-titlepage-recto-style
      font-size: (HSIZE 4)
      line-spacing: (* (HSIZE 4) %line-spacing-factor%)
      space-before: (* (HSIZE 4) %head-before-factor%)
      quadding: %division-subtitle-quadding%
      keep-with-next?: #t
      (process-children-trim)))

;;  (element surname
;;    (make paragraph
;;      use: book-titlepage-recto-style
;;      quadding: %division-title-quadding%
;;      (process-children)))

  (element title 
    (make paragraph
      use: book-titlepage-recto-style
      font-size: (HSIZE 8)
      line-spacing: (* (HSIZE 8) %line-spacing-factor%)
      space-before: (* (HSIZE 8) %head-before-factor%)
      space-after: (* (HSIZE 2) %line-spacing-factor%)
      quadding: %division-title-quadding%
      keep-with-next?: #t
      heading-level: (if %generate-heading-level% 1 0)
      (with-mode title-mode
	(process-children-trim))))

)



