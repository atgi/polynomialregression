<?php

namespace ATGI\PolynomialRegression;

class PolynomialRegression {
	/** @internal Array of running sum of x^n. */
	private $xPowers;

	/** @internal Array of running sum of x*y powers. */
	private $xyPowers;

	/** @internal Number of coefficients. */
	private $numberOfCoefficient;

	/** @internal Array of forcing terms. */
	private $forcedValue;

	/** @internal The index of the current element.  Basiclly a count of data added. */
	private $index = 0;

	/** @internal The weighting interface.  NULL is unused. */
	private $weightingInterface = null;

	/**
	 * Constructor
	 *
	 * Create new class.
	 *
	 * @param $numberOfCoefficient Number of coefficients in polynomial (degree
	 *   of polynomial + 1).
	 */
	public function __construct( $numberOfCoefficient = 3 ) {
		$this->numberOfCoefficient = $numberOfCoefficient;
		$this->reset();

	} // __construct

	/**
	 * Reset data.
	 *
	 * Clear all internal data and prepare for new calculation.
	 * Must be called *after* setNumberOfCoefficient if number of coefficients has
	 * changed.
	 */
	public function reset() {
		$this->forcedValue        = array();
		$this->xPowers            = array();
		$this->xyPowers           = array();
		$this->index              = 0;
		$this->weightingInterface = null;

		$squares = ( $this->numberOfCoefficient - 1 ) * 2;

		// Initialize power arrays.
		for ( $index = 0; $index <= $squares; ++ $index ) {
			$this->xPowers[ $index ]  = 0;
			$this->xyPowers[ $index ] = 0;
		}

	} // reset

	/**
	 * Set degree (Deprecated).
	 *
	 * This is the maximum number of coefficients the polynomial function that
	 * will be calculate.  Note that the request for coefficients can be lower
	 * then this value.  If number is higher, data must be reset and
	 * added again.
	 *
	 * @param int $numberOfCoefficient Number of coefficients.
	 *
	 * @deprecated Deprecated in version 1.1 because the name doesn't reflect what
	 *   the function actually does.  Use 'setNumberOfCoefficient' instead--operations
	 *   are identical.
	 */
	public function setDegree( $numberOfCoefficient ) {
		$this->numberOfCoefficient = $numberOfCoefficient;

	} // setDegree

	/**
	 * Set the number of coefficients to calculate.
	 *
	 * This is the maximum number of coefficients the polynomial function that
	 * will be calculate.  Note that the request for coefficients can be lower
	 * then this value.  If number is higher, data must be reset and
	 * added again.
	 *
	 * @param int $numberOfCoefficient Number of coefficients.
	 *
	 * @since Version 1.1
	 */
	public function setNumberOfCoefficient( $numberOfCoefficient ) {
		$this->numberOfCoefficient = $numberOfCoefficient;

	} // setNumberOfCoefficient

	/**
	 * Get the number of coefficients to calculate.
	 *
	 * Returns the number of coefficients calculated.
	 * @return int Number of coefficients.
	 * @since Version 1.1
	 */
	public function getNumberOfCoefficient( $numberOfCoefficient ) {
		return $this->numberOfCoefficient;

	} // getnumberOfCoefficient

	/**
	 * Set a forced coefficient.
	 *
	 * Force a coefficient to be assumed a specific value for the calculation.
	 * Most often used to force an offset of zero, but can be used to set any
	 * known coefficient for a set of data.
	 *
	 * @param int $coefficient Which coefficient to force.
	 * @param float $value Value to force this coefficient.
	 *
	 * @since Version 1.0
	 */
	public function setForcedCoefficient( $coefficient, $value ) {
		$this->forcedValue[ $coefficient ] = $value;

	} // setForcedCoefficient

	/**
	 * Get a forced coefficient.
	 *
	 * Get a previously set forced coefficient.
	 *
	 * @param int $coefficient Which coefficient.
	 *
	 * @return float Value of this force this coefficient.  Null if the
	 *   coefficient isn't being forced.
	 * @since Version 1.1
	 */
	public function getForcedCoefficient( $coefficient, $value ) {
		$result = null;
		if ( isset( $this->forcedValue[ $coefficient ] ) ) {
			$result = $this->forcedValue[ $coefficient ];
		}

		return $result;

	} // getForcedCoefficient

	/**
	 * Set a weighting interface.
	 *
	 * The regression can be weighted on a per-index basis using a weighting
	 * interface.  An instance of this interface can be set here.
	 *
	 * @param WeightingInterface $weightingInterface Instance of weighting system
	 *   to be used.
	 *
	 * @since Version 1.2
	 */
	public function setWeighting( WeightingInterface $weightingInterface ) {
		$this->weightingInterface = $weightingInterface;
	}

	/**
	 * Get the weighting interface.
	 *
	 * Return the current weighting interface being used.  Returns NULL if no
	 * interface is used.
	 *
	 * @param WeightingInterface $weightingInterface Instance of weighting system
	 *   to be used.
	 *
	 * @since Version 1.2
	 */
	public function getWeighting() {
		return $this->weightingInterface;
	}

	/**
	 * Add data.
	 *
	 * Add a data point to calculation.
	 *
	 * @param float $x Some real value.
	 * @param float $y Some real value corresponding to $x.
	 */
	public function addData( $x, $y ) {
		$squares = ( $this->numberOfCoefficient - 1 ) * 2;

		// Get weighting term for this index.
		$this->index += 1;
		$weight = null;
		if ( null !== $this->weightingInterface ) {
			$weight = $this->weightingInterface->getWeight( $this->index );
		}

		// Remove the effect of the forced coefficient from this value.
		foreach ( $this->forcedValue as $coefficient => $value ) {
			$sub = bcpow( $x, $coefficient );
			$sub = bcmul( $sub, $value );
			$y   = bcsub( $y, $sub );
		}

		// Accumulation of $x raised to the power of the loop iteration.
		// $xSum = pow( $x, $index ) starting with pow( $x, 0 ) which is 1.
		$xSum = 1;

		// Accumulate new data to power sums.
		for ( $index = 0; $index <= $squares; ++ $index ) {
			$accumulator = $xSum;

			// Add weighting term (if applicable).
			if ( null !== $weight ) {
				$accumulator = bcmul( $accumulator, $weight );
			}

			$this->xPowers[ $index ] =
				bcadd( $this->xPowers[ $index ], $accumulator );

			$this->xyPowers[ $index ] =
				bcadd( $this->xyPowers[ $index ], bcmul( $y, $accumulator ) );

			$xSum = bcmul( $xSum, $x );
		}

	} // addData

	/**
	 * Get coefficients.
	 *
	 * Calculate and return coefficients based on current data.
	 *
	 * @param int $numberOfCoefficient Integer value of the degree polynomial desired.  Default
	 *    is -1 which is the max number of coefficients set by class.
	 *
	 * @return array Array of coefficients (as BC strings).
	 */
	public function getCoefficients( $numberOfCoefficient = - 1 ) {
		// If no number of coefficients specified, use standard.
		if ( $numberOfCoefficient == - 1 ) {
			$numberOfCoefficient = $this->numberOfCoefficient;
		}

		// Build a matrix.
		// The matrix is made up of the sum of powers.  So if the number represents the power,
		// the matrix will look like this for a 4th degree polynomial:
		//     [ 0 1 2 3 4 ]
		//     [ 1 2 3 4 5 ]
		//     [ 2 3 4 5 6 ]
		//     [ 3 4 5 6 7 ]
		//     [ 4 5 6 7 8 ]
		//
		$matrix = array();
		for ( $row = 0; $row < $numberOfCoefficient; ++ $row ) {
			$matrix[ $row ] = array();
			for ( $column = 0; $column < $numberOfCoefficient; ++ $column ) {
				$matrix[ $row ][ $column ] =
					$this->xPowers[ $row + $column ];
			}
		}

		// Create augmented matrix by adding X*Y powers.
		for ( $row = 0; $row < $numberOfCoefficient; ++ $row ) {
			$matrix[ $row ][ $numberOfCoefficient ] = $this->xyPowers[ $row ];
		}

		// Add in the forced coefficients.  This is done by nulling the row and column
		// for each forced coefficient.  For example, a 3th degree polynomial
		// matrix with have the 2nd coefficient set to F:
		//       [ a b c d w ]      [ a 0 c d w ]
		//       [ b c d e x ]  ->  [ 0 1 0 0 F ]
		//       [ c d e f y ]      [ c 0 e f y ]
		//       [ d e f g z ]      [ d 0 f g z ]
		foreach ( $this->forcedValue as $coefficient => $value ) {
			for ( $index = 0; $index < $numberOfCoefficient; ++ $index ) {
				$matrix[ $index ][ $coefficient ] = "0";
				$matrix[ $coefficient ][ $index ] = "0";
			}

			$matrix[ $coefficient ][ $coefficient ]         = "1";
			$matrix[ $coefficient ][ $numberOfCoefficient ] = $value;
		}

		// Determine number of rows in matrix.
		$rows = count( $matrix );

		// Initialize done.
		$isDone = array();
		for ( $column = 0; $column < $rows; ++ $column ) {
			$isDone[ $column ] = false;
		}

		// This loop will result in an upper-triangle matrix with the
		// diagonals all 1--the first part of row-reduction--using 2
		// elementary row operations: multiplying a row by a scalar, and
		// subtracting a row by a multiple of an other row.
		// NOTE: This loop can be done out-of-order.  That is, the first
		// row may not begin with the first term.  Order is tracked in the
		// "order" array.
		$order = array();
		for ( $column = 0; $column < $rows; ++ $column ) {
			// Find a row to work with.
			// A row that has a term in this column, and has not yet been
			// reduced.
			$activeRow = 0;
			while ( ( ( 0 == $matrix[ $activeRow ][ $column ] )
			          || ( $isDone[ $activeRow ] ) )
			        && ( $activeRow < $rows ) ) {
				++ $activeRow;
			}

			// Do we have a term in this row?
			if ( $activeRow < $rows ) {
				// Remember the order.
				$order[ $column ] = $activeRow;

				// Normalize row--results in the first term being 1.
				$firstTerm = $matrix[ $activeRow ][ $column ];
				for ( $subColumn = $column; $subColumn <= $rows; ++ $subColumn ) {
					$matrix[ $activeRow ][ $subColumn ] =
						bcdiv( $matrix[ $activeRow ][ $subColumn ], $firstTerm );
				}

				// This row is finished.
				$isDone[ $activeRow ] = true;

				// Subtract the active row from all rows that are not finished.
				for ( $row = 0; $row < $rows; ++ $row ) {
					if ( ( ! $isDone[ $row ] )
					     && ( 0 != $matrix[ $row ][ $column ] )
					) {
						// Get first term in row.
						$firstTerm = $matrix[ $row ][ $column ];
						for ( $subColumn = $column; $subColumn <= $rows; ++ $subColumn ) {
							$accumulator                  = bcmul( $firstTerm, $matrix[ $activeRow ][ $subColumn ] );
							$matrix[ $row ][ $subColumn ] =
								bcsub( $matrix[ $row ][ $subColumn ], $accumulator );
						}
					}
				}
			}
		}

		// Reset done.
		for ( $row = 0; $row < $rows; ++ $row ) {
			$isDone[ $row ] = false;
		}

		$coefficients = array();

		// Back-substitution.
		// This will solve the matrix completely, resulting in the identity
		// matrix in the x-locations, and the coefficients in the last column.
		//   | 1  0  0 ... 0  c0 |
		//   | 0  1  0 ... 0  c1 |
		//   | .  .  .     .   . |
		//   | .  .  .     .   . |
		//   | 0  0  0 ... 1  cn |
		for ( $column = ( $rows - 1 ); $column >= 0; -- $column ) {
			// The active row is based on order.
			$activeRow = $order[ $column ];

			// The active row is now finished.
			$isDone[ $activeRow ] = true;

			// For all rows not finished...
			for ( $row = 0; $row < $rows; ++ $row ) {
				if ( ! $isDone[ $row ] ) {
					$firstTerm = $matrix[ $row ][ $column ];

					// Back substitution.
					for ( $subColumn = $column; $subColumn <= $rows; ++ $subColumn ) {
						$accumulator                  =
							bcmul( $firstTerm, $matrix[ $activeRow ][ $subColumn ] );
						$matrix[ $row ][ $subColumn ] =
							bcsub( $matrix[ $row ][ $subColumn ], $accumulator );
					}
				}
			}

			// Save this coefficient for the return.
			$coefficients[ $column ] = $matrix[ $activeRow ][ $rows ];
		}

		// Coefficients are stored backward, so sort them.
		ksort( $coefficients );

		// Return the coefficients.
		return $coefficients;

	} // getCoefficients

	/**
	 * Interpolate
	 *
	 * Return y point for given x and coefficient set.  Function is static as it
	 * does not require any instance data to operate.
	 *
	 * @param array $coefficients Coefficients as calculated by 'getCoefficients'.
	 * @param float $x X-coordinate from which to calculate Y.
	 *
	 * @return float Y-coordinate (as floating-point).
	 */
	static public function interpolate( $coefficients, $x ) {
		$numberOfCoefficient = count( $coefficients );

		$y = 0;
		for ( $coefficentIndex = 0; $coefficentIndex < $numberOfCoefficient; ++ $coefficentIndex ) {
			// y += coefficients[ coefficentIndex ] * x^coefficentIndex
			$y =
				bcadd
				(
					$y,
					bcmul
					(
						$coefficients[ $coefficentIndex ],
						bcpow( $x, $coefficentIndex )
					)
				);
		}

		return floatval( $y );

	} // interpolate

} // Class
