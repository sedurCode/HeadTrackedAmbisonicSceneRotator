/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
using namespace juce;
#include "AmbisonicRotation.h"
#include "../../ht-api-juce/supperware/HeadMatrix.h"

//==============================================================================
/**
*/
class HeadTrackedAmbisonicSceneRotatorAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    HeadTrackedAmbisonicSceneRotatorAudioProcessor();
    ~HeadTrackedAmbisonicSceneRotatorAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    //==============================================================================
	void setGain(float newGain) { fTargetGain = newGain; };
	void setYaw(float newYaw)   { fTargetYaw = newYaw; };
	void setPitch(float newPitch) { fTargetPitch = newPitch; };
	void setRoll(float newRoll) { fTargetRoll = newRoll; };
    void setRotationMatrix(const HeadMatrix& headMatrix) { ambisonicRotator.updateRotation(headMatrix); };

private:
    //==============================================================================

	AmbisonicRotation ambisonicRotator;

	// Targets
	float fTargetGain = 1.f;
	float fTargetYaw = 0.f;
	float fTargetPitch = 0.f;
	float fTargetRoll = 0.f;
	// Current Values
	float fCurrentGain = 0.f;
	float fCurrentYaw = 0.f;
	float fCurrentPitch = 0.f;
	float fCurrentRoll = 0.f;
	// Steps
	float fGainStep = 0.f;
	float fYawStep = 0.f;
	float fPitchStep = 0.f;
	float fRollStep = 0.f;
	// Constants
	const float fMaxDeltaStep = 1.f / 480.f;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (HeadTrackedAmbisonicSceneRotatorAudioProcessor)
};
