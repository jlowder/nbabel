;; This version uses high-speed vector functions based on arrays. Acceleration values are cached as
;; members of the star class for efficient implementation of the leapfrog method.

(declaim (optimize (speed 3) (debug 0) (safety 0))
         (inline make-3v v+v v-v vlen^2 vlen s*v force 3v/0 3v/1 3v/2))

(setf *read-default-float-format* 'double-float)

(defconstant +dt+ 1d-3) ; global timestep

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; a 3-vector type
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(deftype 3v () '(simple-array double-float (3)))

(defun make-3v (x y z)
  (make-array 3 :element-type 'double-float :initial-contents (list x y z)))

(defun 3v/0 (v)
  (declare (type 3v v))
  (aref v 0))

(defun 3v/1 (v)
  (declare (type 3v v))
  (aref v 1))

(defun 3v/2 (v)
  (declare (type 3v v))
  (aref v 2))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; some vector utility functions based on the 3v type
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun v+v (a b)
  "add two vectors"
  (declare (type 3v a b))
  (make-3v (+ (3v/0 a) (3v/0 b))
           (+ (3v/1 a) (3v/1 b))
           (+ (3v/2 a) (3v/2 b))))

(defun v-v (a b)
  "subtract two vectors"
  (declare (type 3v a b))
  (make-3v (- (3v/0 a) (3v/0 b))
           (- (3v/1 a) (3v/1 b))
           (- (3v/2 a) (3v/2 b))))

(defun vlen^2 (v)
  "length of a vector, squared"
  (declare (type 3v v))
  (+ (* (3v/0 v) (3v/0 v))
     (* (3v/1 v) (3v/1 v))
     (* (3v/2 v) (3v/2 v))))

(defun vlen (v)
  "length of a vector"
  (declare (type 3v v))
  (sqrt (vlen^2 v)))

(defun s*v (s v)
  "scalar times vector"
  (declare (type 3v v))
  (make-3v (* s (3v/0 v))
           (* s (3v/1 v))
           (* s (3v/2 v))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; a class representing a star (i.e. particle)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass star ()
  ((pos :type 3v :accessor pos :initarg :pos :initform (make-3v 0d0 0d0 0d0))
   (vel :type 3v :accessor vel :initarg :vel :initform (make-3v 0d0 0d0 0d0))
   (mass :accessor mass :initarg :mass :initform 0d0)
   (a0 :type 3v :accessor a0 :initarg :a0 :initform (make-3v 0d0 0d0 0d0))
   (acc :type 3v :accessor acc :initarg :acc :initform (make-3v 0d0 0d0 0d0))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; a class representing a star cluster
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass cluster ()
  ((stars :accessor stars :initarg :stars :initform '())))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; generic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric potential-energy (c))

(defgeneric kinetic-energy (c))

(defgeneric energy (c))

(defgeneric update-velocities (c dt))

(defgeneric update-positions (c dt))

(defgeneric update-accelerations (c))

(defgeneric load-stdin (c))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; specialized methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod potential-energy ((c cluster))
  "potential energy calculation"
  (let* ((stars (stars c))
         (n (length stars)))
    (loop for i below n sum
         (loop for j from (1+ i) below n
            as star-i = (elt stars i)
            as star-j = (elt stars j)
            sum (/ (* -1 (mass star-i) (mass star-j)) (vlen (v-v (pos star-j) (pos star-i))))))))

(defmethod kinetic-energy ((c cluster))
  "kinetic energy calculation"
  (loop for star in (stars c)
     sum (* 1/2 (mass star) (vlen^2 (vel star)))))

(defmethod energy ((c cluster))
  "return total energy, potential energy, and kinetic energy"
  (let ((p (potential-energy c))
        (k (kinetic-energy c)))
    (values (+ p k) p k)))

(defmethod update-positions ((c cluster) dt)
  (loop for star in (stars c)
     as vel = (vel star)
     as acc = (acc star)
     do
       (setf (a0 star) acc)
       (setf (pos star) (v+v (pos star) (v+v (s*v dt vel) (s*v (* 1/2 dt dt) acc))))))

(defmethod update-velocities ((c cluster) dt)
  (loop for star in (stars c)
     as ai = (v+v (a0 star) (acc star))
     do
       (setf (vel star) (v+v (vel star) (s*v (* 1/2 dt) ai)))
       (setf (a0 star) (acc star))))

(defun force (ri rj mj)
  "calculate the force between two particles"
  (declare (type 3v ri rj))
  (let* ((p (v-v ri rj))
         (c (* mj (expt (vlen^2 p) -1.5))))
    (s*v c p)))

(defmethod update-accelerations ((c cluster))
  (loop for star in (stars c) do (setf (acc star) (make-3v 0d0 0d0 0d0)))
  (let* ((stars (stars c))
         (n (length stars)))
    (loop for i below n
       as star-i = (elt stars i)
       do (loop for j below n
             as star-j = (elt stars j)
             do (when (not (= i j))
                  (setf (acc star-i) (v-v (acc star-i) (force (pos star-i) (pos star-j) (mass star-j)))))))))

(defun stdin->strlist ()
  "read a text file as a list of strings"
  (loop for line = (read-line *standard-input* nil)
     while line
     collect line))

(defun split-spaces (s &optional acc)
  (let ((o (search " " s)))
    (if (null o)
        (reverse (remove-if (lambda (x) (string= "" x)) (cons s acc)))
        (split-spaces (string-trim '(#\ ) (subseq s o)) (cons (subseq s 0 o) acc)))))

(defmethod load-stdin ((c cluster))
  "initialize the cluster from stdin"
  (setf (stars c)
        (loop for line in (remove-if (lambda (x) (< (length x) 3))
                                     (mapcar #'split-spaces (stdin->strlist)))
           as parsed = (mapcar (lambda (x) (with-input-from-string (l x) (read l))) line)
           as star = (make-instance 'star
                                    :mass (elt parsed 1)
                                    :pos (apply #'make-3v (subseq parsed 2 5))
                                    :vel (apply #'make-3v (subseq parsed 5)))
           collect star)))

(defmethod print-object ((obj star) out)
  (print-unreadable-object (obj out :type t)
    (format out "(~a ~a ~a) (~a ~a ~a) ~a"
            (3v/0 (pos obj))
            (3v/1 (pos obj))
            (3v/2 (pos obj))
            (3v/0 (vel obj))
            (3v/1 (vel obj))
            (3v/2 (vel obj))
            (mass obj))))

(defmethod print-object ((obj cluster) out)
  (print-unreadable-object (obj out :type t)
    (progn
      (format out "cluster: ~%")
      (loop for star in (stars obj) do (format out " ~S~%" star)))))

(defun lineout (time ke pe total err)
  (format t "~A~24T~A~48T~A~72T~A~96T~A~%" time ke pe total err))

(defun summary (time tot ke pe e0)
  (lineout time ke pe tot (/ (- tot e0) e0)))

(defun main ()
  (let ((brk "=======================")
        (cluster (make-instance 'cluster))
        (tend 1d0))
    (load-stdin cluster)
    ;;(format t "~S~%" cluster)
    (update-accelerations cluster)
    (multiple-value-bind (e0 pe ke) (energy cluster)
      (declare (ignore pe ke))
      (lineout "Time" "Kinetic Energy" "Potential Energy" "Total Energy" "Energy Error")
      (lineout brk brk brk brk brk)
      (loop for time from 0d0 to (+ +dt+ tend) by +dt+
         as k = 0 then (1+ k)
         do
           (update-positions cluster +dt+)
           (update-accelerations cluster)
           (update-velocities cluster +dt+)
           (when (zerop (mod k 10))
             (multiple-value-bind (e pe ke) (energy cluster)
               (summary time e ke pe e0)))))))

