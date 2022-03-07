/*
  ==============================================================================

    AmbisonicRotation.h
    Created: 5 Feb 2022 2:40:59pm
    Author:  sedur

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
using namespace juce::dsp;
#include "../../ht-api-juce/supperware/HeadMatrix.h"
class AmbisonicRotation
{
public:
	//void init();
	AmbisonicRotation();
	~AmbisonicRotation();

	void process(juce::AudioSampleBuffer& buffer);
	void calcRotationMatrix(const int order);
	void updateEulerRPY(float r, float p, float y);
	void updateEulerYPR(float y, float p, float r);
	void updateRotation(const HeadMatrix& headMatrix);

	using matrix_float_ptr = juce::dsp::Matrix<float>*;

	float radiansToDegrees(float radians) { return radians * (float(180) / float(3.14159265358979323846264338327950288)); }
	float degreesToRadians(float degrees) { return degrees * (float(3.14159265358979323846264338327950288) / float(180)); }

private:

	float yaw = 0.0f;
	float pitch = 0.0f;
	float roll = 0.0f;
	int cachedOrder = 1;
	bool rotationSequence = true;

	juce::Atomic<bool> rotationParamsHaveChanged{ true };

	juce::AudioSampleBuffer copyBuffer;

	juce::OwnedArray<Matrix<float>> orderMatrices;
	juce::OwnedArray<Matrix<float>> orderMatricesCopy;

	double P(int i, int l, int a, int b, Matrix<float>& R1, Matrix<float>& Rlm1);
	double U(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1);
	double V(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1);
	double W(int l, int m, int n, Matrix<float>& Rone, Matrix<float>& Rlm1);




};