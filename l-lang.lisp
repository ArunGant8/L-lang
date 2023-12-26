;; A tutorial quantum interpreter in 150 lines of lisp

;; We are writing an interpreter for the quantum language L
;; which is, in some sense, the simplest non-trivial quantum language.

;; We implement the semantics operationally, using an abstract machine
;; whose state is an ordered pair consisting of a quantum state and
;; an n-bit measurement register (n > 0, integer) as follows:

(defstruct machine
  quantum-state
  measurement-register)

;; We represent a quantum state for an n-qubit system as an array
;; of 2^n complex numbers. An index i into the array represents a
;; probability amplitude psi_i, which is the scalar component of |i>

;; The following function allocates a new quantum state of n qubits,
;; initialized to be in the |...000> state

(defun make-quantum-state (n)
  (let ((s (make-array (expt 2 n) :initial-element 0.0d0)))
    (setf (aref s 0) 1.0d0)
    s))

;; We can find the number of qubits a state represents, or an operator
;; acts on. It reduces to determining the number of qubits that a
;; dimension represents, equivalent to taking a binary logarithm.
;; In Lisp, we can do this by computing the number of bits an integer
;; takes to represent using integer-length.

(defun dimension-qubits (d)
  (1- (integer-length d)))  ;; we subtract one since the number 2^n
                            ;; is a 1 followed by n 0's, with a total
                            ;; length of n+1

;; Matrix-vector and matrix-matrix multiplication, adapated to
;; (and named for) a quantum computing setting. Same old algorithms,
;; nothing fancy. These are the principal ways in which our quantum
;; program executes (or our quantum state EVOLVES):

(defun apply-operator (matrix column)
  (let* ((matrix-size (array-dimension matrix 0))
	 (result (make-array matrix-size :initial-element 0.0d0)))
    (dotimes (i matrix-size)
      (let ((element 0))
	(dotimes (j matrix-size)
	  (incf element (* (aref matrix i j) (aref column j))))
	(setf (aref result i) element)))
    (replace column result)))

(defun compose-operators (A B)
  (destructuring-bind (m n) (array-dimensions A)
    (let* ((l (array-dimension B 1))
	   (result (make-array (list m l) :initial-element 0)))
      (dotimes (i m result)
	(dotimes (k l)
	  (dotimes (j n)
	    (incf (aref result i k)
		  (* (aref A i j)
		     (aref B j k)))))))))

;; MEASUREMENT: Implemented in two steps: as the sampling of
;; the state followed by its collapse.

(defun observe (machine)
  (let ((b (sample (machine-quantum-state machine))))
    (collapse (machine-quantum-state machine) b)
    (setf (machine-measurement-register machine) b)
    machine))

;; Now we define SAMPLE and COLLAPSE

;; We sample by first sampling from (0, 1),
;; then finding the partial sums till we find
;; an interval where our observation falls.

(defun sample (state)
  (let ((r (random 1.0d0)))
    (dotimes (i (length state))
      (decf r (expt (abs (aref state i)) 2))
      (when (minusp r) (return i)))))

;; Collapsing to |k> is simply zeroing out the
;; array and setting psi_k = 1

(defun collapse (state basis-element)
  (fill state 0.0d0)
  (setf (aref state basis-element) 1.0d0))

;; GATES: unitary matrices

;; Identity gate:
(defparameter +I+ #2A((1 0)
		      (0 1)))

;; A general function apply-gate to apply any kind of gate
;; on any collection of qubits for any quantum state

(defun apply-gate (state U qubits)
  (assert (= (length qubits) (dimension-qubits (array-dimension U 0))))
  (if (= 1 (length qubits))
      (%apply-1Q-gate state U (first qubits))
      (%apply-nQ-gate state U qubits)))

;; LIFTING: A process by which we generate large-dimensional
;; unitary matrices using a one qubit gate matrix, the qubit index,
;; and the size of the machine. This is essentially an application
;; of the Kronecker product to construct a larger operator that
;; works on the tensor product of two smaller spaces, from two
;; operators that operate on those spaces.

;; So, first we define the Kronecker product itself

(defun kronecker-multiply (A B)
  (destructuring-bind (m n) (array-dimensions A)
    (destructuring-bind (p q) (array-dimensions B)
      (let ((result (make-array (list (* m p) (* n q)))))
	(dotimes (i m result)
	  (dotimes (j n)
	    (let ((Aij (aref A i j))
		  (y (* i p))
		  (x (* j q)))
	      (dotimes (u p)
		(dotimes (v q)
		  (setf (aref result (+ y u) (+ x v))
			(* Aij (aref B u v))))))))))))

;; Then we define a way to "exponentiate" a Kronecker product

(defun kronecker-expt (U n)
  (cond
    ((< n 1) #2A((1)))
    ((= n 1) U)
    (t (kronecker-multiply (kronecker-expt U (1- n)) U))))

;; Now we can write LIFT:
;; This takes care of single qubit gates as well as
;; gates on ADJACENT qubits. (The logic is similar)

(defun lift (U i n)
  (let ((left (kronecker-expt +I+ (- n i (dimension-qubits
					  (array-dimension U 0)))))
	(right (kronecker-expt +I+ i)))
    (kronecker-multiply left (kronecker-multiply U right))))

;; To deal with multi-qubit gates with non-adjacent qubits,
;; we first SWAP the non-adjacent qubits till we have them
;; adjacent, then apply the previous method, then perform
;; the sequence of permutations in reverse.

;; So we start by defining a SWAP operator:
(defparameter +SWAP+ #2A((1 0 0 0)
			 (0 0 1 0)
			 (0 1 0 0)
			 (0 0 0 1)))

;; The above operator (and its lifted variants) can swap
;; two index-adjacent qubits. But this is enough: Any permutation
;; can be decomposed into a composition of swaps, and every
;; swap can be decomposed into a series of adjacent transpositions.

;; We define a function that takes a permutation written as a list,
;; and converts it to a list of (possibly non-adjacent) transpositions
;; to be applied L-to-R

(defun permutation-to-transpositions (permutation)
  (let ((swaps nil))
    (dotimes (dest (length permutation) (nreverse swaps))
      (let ((src (elt permutation dest)))
	(loop :while (< src dest) :do
	  (setf src (elt permutation src)))
	(cond
	  ((< src dest) (push (cons src dest) swaps))
	  ((> src dest) (push (cons dest src) swaps)))))))

;; Next we convert these transpositions to adjacent transposition indexes.

(defun transposition-to-adjacent-transpositions (transpositions)
  (flet ((expand-cons (c)
	   (if (= 1 (- (cdr c) (car c)))
	       (list (car c))
	       (let ((trans (loop :for i :from (car c) :below (cdr c)
				  :collect i)))
		 (append trans (reverse (butlast trans)))))))
    (mapcan #'expand-cons transpositions)))

;; We are finally ready to write the function to apply
;; a multi-qubit gate to a multi-qubit state

(defun %apply-nQ-gate (state U qubits)
  (let ((n (dimension-qubits (length state))))
    (labels ((swap (i)
	       (lift +swap+ i n))
	     (transpositions-to-operator (trans)
	       (reduce #'compose-operators trans :key #'swap)))
      (let* ((U01 (lift U 0 n))
	     (from-space (append (reverse qubits)
				 (loop :for i :below n
				       :when (not (member i qubits))
					 :collect i)))
	     (trans (transpositions-to-adjacent-transpositions
		     (permutation-to-transpositions
		      from-space)))
	     (to->from (transpositions-to-operator trans))
	     (from->to (transpositions-to-operator (reverse trans)))
	     (Upq (compose-operators to->from
				     (compose-operators U01
							from->to))))
	(apply-operator Upq state)))))

;; Now for the driver loop:

(defun run-quantum-program (qprog machine)
  (loop :for (instruction . payload) :in qprog
	:do (ecase instruction
	      ((GATE)
	       (destructuring-bind (gate &rest qubits) payload
		 (apply-gate (machine-quantum-state machine) gate qubits)))
	      ((MEASURE)
	       (observe machine)))
	    :finally (return machine)))
