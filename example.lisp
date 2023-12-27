;; We define the BELL STATE:

;; H operator:

(defparameter +H+ (make-array '(2 2) :initial-contents (let ((s (/ (sqrt 2))))
								(list (list s s)
								      (list s (- s))))))

;; CNOT operator:

(defparameter +CNOT+ #2A((1 0 0 0)
			 (0 1 0 0)
			 (0 0 0 1)
			 (0 0 1 0)))

;; Now the Bell State is synthesized as:

(defun bell (p q)
  `((GATE ,+H+ ,p)
    (GATE ,+CNOT+ ,p ,q)))

;; Now the Greenberger-Horne-Zellinger state,
;; a generalization of the Bell state for more
;; than two qubits.

(defun ghz (n)
  (cons `(GATE ,+H+ 0)
	(loop :for q :below (1- n)
	      :collect `(GATE ,+CNOT+ ,q ,(1+ q)))))

;; The quantum Fourier transform:

;; - the CPHASE gate

(defun cphase (angle)
  (make-array '(4 4) :initial-contents `((1 0 0 0)
					 (0 1 0 0)
					 (0 0 1 0)
					 (0 0 0 ,(cis angle)))))

(defun qft (qubits)
  (labels ((bit-reversal (qubits)
			 (let ((n (length qubits)))
			   (if (< n 2)
			       nil
			     (loop :repeat (floor n 2)
				   :for qs :in qubits
				   :for qe :in (reverse qubits)
				   :collect `(GATE ,+swap+ ,qs ,qe)))))
	   (%qft (qubits)
		 (destructuring-bind (q . qs) qubits
				     (if (null qs)
					 (list `(GATE ,+H+ ,q))
				       (let ((cR (loop :with n := (1+ (length qs))
						       :for i :from 1
						       :for qi :in qs
						       :for angle := (/ pi (expt 2 (- n i)))
						       :collect `(GATE ,(cphase angle) ,q ,qi))))
					 (append
					  (qft qs)
					  cR
					  (list `(GATE ,+H+ ,q))))))))
	  (append (%qft qubits) (bit-reversal qubits))))
