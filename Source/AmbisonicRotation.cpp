/*
  ==============================================================================

    AmbisonicRotation.cpp
    Created: 5 Feb 2022 2:40:59pm
    Author:  sedur

  ==============================================================================
  This code implements rotation stategies in higher ambisonics taken from the paper:
  Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion
  by Joseph Ivanic and Klaus Ruedenberg (1996)
  With the corrections and additions proposed in 1998
  Published in The Journal of Phyisical Chemistry. 
  The application of this approach was adapted from Archonits Politis HoA Matlab Library (2015)
  https://github.com/polarch/Higher-Order-Ambisonics
  Which also relies on the spherical harmonic transform and VBAP libraries from the same author
  https://github.com/polarch/Spherical-Harmonic-Transform
  https://github.com/polarch/Vector-Base-Amplitude-Panning

*/

#include "AmbisonicRotation.h"

AmbisonicRotation::AmbisonicRotation()
{
	orderMatrices.add(new Matrix<float>(0, 0)); // 0th
	orderMatricesCopy.add(new Matrix<float>(0, 0)); // 0th

	for (int l = 1; l <= 7; ++l)
	{
		const int nCh = (2 * l + 1);
		auto elem = orderMatrices.add(new Matrix<float>(nCh, nCh));
		elem->clear();
		auto elemCopy = orderMatricesCopy.add(new Matrix<float>(nCh, nCh));
		elemCopy->clear();
	}
}

AmbisonicRotation::~AmbisonicRotation()
{}

void AmbisonicRotation::process(juce::AudioSampleBuffer& buffer)
{
	// Dynamically set up based on number of channels
	// TODO: Make this switchable with a GUI element instead
	const int numberOfChannels = buffer.getNumChannels();
	const int actualOrder = floor(sqrt(numberOfChannels)) - 1;
	cachedOrder = actualOrder;
	const int bufferLength = buffer.getNumSamples();
	bool newRotationMatrix = false;

	// Perform new rotation if needed
	if (rotationParamsHaveChanged.get())
	{
		newRotationMatrix = true;
		calcRotationMatrix(actualOrder);
	}

	// Clear buffer working copy
	copyBuffer.clear();
	copyBuffer.setSize(numberOfChannels, bufferLength);

	// make copy of the input buffer
	for (int ch = 0; ch < numberOfChannels; ++ch)
		copyBuffer.copyFrom(ch, 0, buffer, ch, 0, bufferLength);

	// clear all channels except W (mono)
	for (int ch = 1; ch < buffer.getNumChannels(); ++ch)
		buffer.clear(ch, 0, bufferLength);

	// Operating over orders 1 to n
	for (int l = 1; l <= actualOrder; ++l)
	{
		const int offset = l * l;
		const int nCh = 2 * l + 1;

		//Get the rotation matrices for this order
		matrix_float_ptr R = orderMatrices[l];
		matrix_float_ptr Rcopy = orderMatricesCopy[l];

		// over the number of channels in this ambisonic order
		for (int o = 0; o < nCh; ++o)
		{
			// from offset to offset + the number of channels in this order
			const int chOut = offset + o;
			for (int p = 0; p < nCh; ++p)
			{
				// TODO: If we don't have a new rotation matrix, we can do all these operations without costly ramping
				// Sum into this channel of the buffer the content of 
				buffer.addFromWithRamp(chOut,                                 // sum in to the target channel
									   0,                                     // starting at sample 0
									   copyBuffer.getReadPointer(offset + p), // the (offset + p)th channel of content
									   bufferLength,                          // over the buffer length
								       Rcopy->operator() (o, p),              // Ramping from the current position in the target rotation matrix
									   R->operator()     (o, p)               // to the desired rotation of the matrix
				);    
			}
		}
	}

	// make copies for fading between old and new matrices
	if (newRotationMatrix)
	{
		for (int l = 1; l <= actualOrder; ++l) // Across the order
		{
			*orderMatricesCopy[l] = *orderMatrices[l]; // copy the rotation matrices across to update the current state now that we have ramped
		}
	}
}

double AmbisonicRotation::P(int i, int l, int a, int b, Matrix<float>& R1, Matrix<float>& Rlm1)
{
	double ri1  = R1(i + 1, 2);
	double rim1 = R1(i + 1, 0);
	double ri0  = R1(i + 1, 1);

	if (b == -l)
		return ri1 * Rlm1(a + l - 1, 0) + rim1 * Rlm1(a + l - 1, 2 * l - 2);
	else if (b == l)
		return ri1 * Rlm1(a + l - 1, 2 * l - 2) - rim1 * Rlm1(a + l - 1, 0);
	else
		return ri0 * Rlm1(a + l - 1, b + l - 1);
};

double AmbisonicRotation::U(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1)
{
	return P(0, l, m, n, Rone, Rlm1);
}

double AmbisonicRotation::V(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1)
{
	if (m == 0) // If the mth (major) harmonic is zero
	{
		auto p0 = P( 1, l,  1, n, Rone, Rlm1);
		auto p1 = P(-1, l, -1, n, Rone, Rlm1);
		return p0 + p1;
	}
	else if (m > 0) // If the mth (major) harmonic is positive and non-zero
	{
		auto p0 = P(1, l, m - 1, n, Rone, Rlm1);
		if (m == 1) // d = 1;
			return p0 * sqrt(2);
		else // d = 0;
			return p0 - P(-1, l, 1 - m, n, Rone, Rlm1);
	}
	else // If the mth (major) harmonic is negative and non-zero
	{
		auto p1 = P(-1, l, -m - 1, n, Rone, Rlm1);
		if (m == -1) // d = 1;
			return p1 * sqrt(2);
		else // d = 0;
			return p1 + P(1, l, m + 1, n, Rone, Rlm1);
	}
}

/**
*/
double AmbisonicRotation::W(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1)
{
    // If the mth (major) harmonic is positive and non-zero
	if (m > 0)
	{
		auto p0 = P( 1 , l,  m + 1, n, Rone, Rlm1);
		auto p1 = P(-1,  l, -m - 1, n, Rone, Rlm1);
		return p0 + p1;
	}
	else if (m < 0)
	{
		auto p0 = P( 1, l, m - 1, n, Rone, Rlm1);
		auto p1 = P(-1, l, 1 - m, n, Rone, Rlm1);
		return p0 - p1;
	}

	return 0.0;
}

/**
 * 
 */
void AmbisonicRotation::calcRotationMatrix(const int order)
{
	const float yawRadians   = degreesToRadians(yaw) * -1; // TODO: We should probably do inversions evenly for everything
	const float pitchRadians = degreesToRadians(pitch);
	const float rollRadians  = degreesToRadians(roll);

	// Calculate the sines and cosines for each of the rotation elements
	float ca = std::cos(yawRadians);
	float cb = std::cos(pitchRadians);
	float cy = std::cos(rollRadians);
	float sa = std::sin(yawRadians);
	float sb = std::sin(pitchRadians);
	float sy = std::sin(rollRadians);

	// Declare a new rotatation matrix
	Matrix<float> rotMat(3, 3);


	if (rotationSequence == true) // roll -> pitch -> yaw (extrinsic rotations)
	{
		// General rotation matrix factors, see https://en.wikipedia.org/wiki/Rotation_matrix
		rotMat(0, 0) = ca * cb;
		rotMat(1, 0) = sa * cb;
		rotMat(2, 0) = -sb;

		rotMat(0, 1) = ca * sb * sy - sa * cy;
		rotMat(1, 1) = sa * sb * sy + ca * cy;
		rotMat(2, 1) = cb * sy;

		rotMat(0, 2) = ca * sb * cy + sa * sy;
		rotMat(1, 2) = sa * sb * cy - ca * sy;
		rotMat(2, 2) = cb * cy;
	}
	else // yaw -> pitch -> roll (extrinsic rotations)
	{
		rotMat(0, 0) = ca * cb;
		rotMat(1, 0) = sa * cy + ca * sb * sy;
		rotMat(2, 0) = sa * sy - ca * sb * cy;

		rotMat(0, 1) = -sa * cb;
		rotMat(1, 1) = ca * cy - sa * sb * sy;
		rotMat(2, 1) = ca * sy + sa * sb * cy;

		rotMat(0, 2) = sb;
		rotMat(1, 2) = -cb * sy;
		rotMat(2, 2) = cb * cy;
	}

	{
		// Place (and re-order) the rotation matrics in the set of order matrices
		matrix_float_ptr Rl = orderMatrices[1];
		Rl->operator() (0, 0) = rotMat(1, 1);
		Rl->operator() (0, 1) = rotMat(1, 2);
		Rl->operator() (0, 2) = rotMat(1, 0);
		Rl->operator() (1, 0) = rotMat(2, 1);
		Rl->operator() (1, 1) = rotMat(2, 2);
		Rl->operator() (1, 2) = rotMat(2, 0);
		Rl->operator() (2, 0) = rotMat(0, 1);
		Rl->operator() (2, 1) = rotMat(0, 2);
		Rl->operator() (2, 2) = rotMat(0, 0);

	}


	// Calculate the higher order rotation matrices
	// This process takes advantage of the recurrent 
	// Operating from second order upwards to the target order (as we already have the first order matrix)
	for (int l = 2; l <= order; ++l)
	{
	    // Take the basis matrix for zero order
		matrix_float_ptr Rone = orderMatrices[1];
		// Take the matrix of the previous order
		matrix_float_ptr Rlm1 = orderMatrices[l - 1];
		// Get a pointer to the matrix we want to fulfill
		matrix_float_ptr Rl   = orderMatrices[l];
		// Working across the square of the range of spherical harmonics in the order
		for (int m = -l; m <= l; ++m)
		{
			for (int n = -l; n <= l; ++n)
			{
				const int d = (m == 0) ? 1 : 0; //if we are operating on the 0th harmonic of the order
				double denom;
				// If we are operating on the highest harmonic in the order in the minor loop
				if (abs(n) == l)
			    {
			        denom = (2 * l) * (2 * l - 1); // We don't want a divide by zero! At n=l force the denom to be a large number 
			    } else {
			        denom = l * l - n * n;         // Order^2 - minor_harmonic^2
			    }
					
                // These expressions are all driven to 0 (1 - 1) when d is 1 i.e. when operating on the 0th harmonic
                //
                // sqrt of the quotient of the distance between the order and the major and minor spherical harmonics
				double u = sqrt((l * l - m * m) / denom); // This number will be relatively small or zero for m=n=
				// 
				double v = sqrt((1.0 + d) * (l + abs(m) - 1.0) * (l + abs(m)) / denom) * (1.0 - 2.0 * d) * 0.5;   // 
				// 
				double w = sqrt((l - abs(m) - 1.0) * (l - abs(m)) / denom) * (1.0 - d) * (-0.5);                  // Driven to 0 when operating on the 0th harmonic
                
				//
				if (u != 0.0)

					u *= U(l, m, n, *Rone, *Rlm1);
				if (v != 0.0)

					v *= V(l, m, n, *Rone, *Rlm1);
				if (w != 0.0)

					w *= W(l, m, n, *Rone, *Rlm1);
				// load the target order matrix with the gain
				Rl->operator() (m + l, n + l) = u + v + w; // 
			}
		}
	}

	rotationParamsHaveChanged = false;
}

void AmbisonicRotation::updateEulerRPY(float r, float p, float y)
{
	if (roll != r || pitch != p || yaw != y)
	{
		rotationSequence = true;
		roll = r;
		pitch = p;
		yaw = y;
		rotationParamsHaveChanged = true;
	}
}

void AmbisonicRotation::updateEulerYPR(float y, float p, float r)
{
	if (yaw != y || pitch != p || roll != r)
	{
		rotationSequence = false;
		yaw = y;
		pitch = p;
		roll = r;
		rotationParamsHaveChanged = true;
	}
}

//void AmbisonicRotation::updateRotation(const HeadMatrix& headMatrix)
//{
//
//	{   // We already have a rotation matrix provided by supperware
//		// Place (and re-order) the rotation matrics in the set of order matrices
//		matrix_float_ptr Rl = orderMatrices[1]; //
//		headMatrix.getMatrix(Rl->operator() (0, 0), 4);  //Rl->operator() (0, 0) = rotMat(1, 1);
//		headMatrix.getMatrix(Rl->operator() (0, 1), 5);   //Rl->operator() (0, 1) = rotMat(1, 2); 
//		headMatrix.getMatrix(Rl->operator() (0, 2), 3);   //Rl->operator() (0, 2) = rotMat(1, 0);
//		headMatrix.getMatrix(Rl->operator() (1, 0), 7);   //Rl->operator() (1, 0) = rotMat(2, 1);
//		headMatrix.getMatrix(Rl->operator() (1, 1), 9);   //Rl->operator() (1, 1) = rotMat(2, 2);
//		headMatrix.getMatrix(Rl->operator() (1, 2), 6);   //Rl->operator() (1, 2) = rotMat(2, 0);
//		headMatrix.getMatrix(Rl->operator() (2, 0), 1);   //Rl->operator() (2, 0) = rotMat(0, 1);
//		headMatrix.getMatrix(Rl->operator() (2, 1), 2);   //Rl->operator() (2, 1) = rotMat(0, 2); 
//		headMatrix.getMatrix(Rl->operator() (2, 2), 0);   //Rl->operator() (2, 2) = rotMat(0, 0); 
//
//	}//headMatrix
//	int order = cachedOrder;
//
//	// Calculate the higher order rotation matrices
//	// This process takes advantage of the recurrent 
//	// Operating from second order upwards to the target order (as we already have the first order matrix)
//	for (int l = 2; l <= order; ++l)
//	{
//		// Take the basis matrix for zero order
//		matrix_float_ptr Rone = orderMatrices[1];
//		// Take the matrix of the previous order
//		matrix_float_ptr Rlm1 = orderMatrices[l - 1];
//		// Get a pointer to the matrix we want to fulfill
//		matrix_float_ptr Rl = orderMatrices[l];
//		// Working across the square of the range of spherical harmonics in the order
//		for (int m = -l; m <= l; ++m)
//		{
//			for (int n = -l; n <= l; ++n)
//			{
//				const int d = (m == 0) ? 1 : 0; //if we are operating on the 0th harmonic of the order
//				double denom;
//				// If we are operating on the highest harmonic in the order in the minor loop
//				if (abs(n) == l)
//				{
//					denom = (2 * l) * (2 * l - 1); // We don't want a divide by zero! At n=l force the denom to be a large number 
//				}
//				else {
//					denom = l * l - n * n;         // Order^2 - minor_harmonic^2
//				}
//
//				// These expressions are all driven to 0 (1 - 1) when d is 1 i.e. when operating on the 0th harmonic
//				//
//				// sqrt of the quotient of the distance between the order and the major and minor spherical harmonics
//				double u = sqrt((l * l - m * m) / denom); // This number will be relatively small or zero for m=n=
//				// 
//				double v = sqrt((1.0 + d) * (l + abs(m) - 1.0) * (l + abs(m)) / denom) * (1.0 - 2.0 * d) * 0.5;   // 
//				// 
//				double w = sqrt((l - abs(m) - 1.0) * (l - abs(m)) / denom) * (1.0 - d) * (-0.5);                  // Driven to 0 when operating on the 0th harmonic
//
//				//
//				if (u != 0.0)
//
//					u *= U(l, m, n, *Rone, *Rlm1);
//				if (v != 0.0)
//
//					v *= V(l, m, n, *Rone, *Rlm1);
//				if (w != 0.0)
//
//					w *= W(l, m, n, *Rone, *Rlm1);
//				// load the target order matrix with the gain
//				Rl->operator() (m + l, n + l) = u + v + w; // 
//			}
//		}
//	}
//
//	rotationParamsHaveChanged = false;
//}