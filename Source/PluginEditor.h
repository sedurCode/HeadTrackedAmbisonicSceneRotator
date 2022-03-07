/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include <JuceHeader.h>
using namespace juce;
#include "PluginProcessor.h"
#include "../../ht-api-juce/supperware/HeadMatrix.h"
#include "../../ht-api-juce/supperware/Tracker.h"
#include "../../ht-api-juce/supperware/midi/midi.h"
#include "../../ht-api-juce/supperware/configpanel/configPanel.h"
#include "../../ht-api-juce/supperware/headPanel/headPanel.h"

//==============================================================================
/**
*/
class HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor  : public juce::AudioProcessorEditor
                                                            , private juce::Slider::Listener
                                                            , public HeadPanel::HeadPanel::Listener
{
public:
    HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor (HeadTrackedAmbisonicSceneRotatorAudioProcessor&);
    ~HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor() override;

    //==============================================================================
    void paint (juce::Graphics&) override;
    void resized() override;
    void trackerChanged(const HeadMatrix& headMatrix) override;
private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    HeadTrackedAmbisonicSceneRotatorAudioProcessor& audioProcessor;
	juce::Slider gainSlider;
	juce::Slider yawSlider;
	juce::Slider pitchSlider;
	juce::Slider rollSlider;
    float intermediateMatrix[9];
    HeadPanel::HeadPanel headPanel;

    Vector3D<float> V;

    void sliderValueChanged(juce::Slider* slider) override;

	float radiansToDegrees(float radians) { return radians * (float(180) / float(3.14159265358979323846264338327950288)); }
	float degreesToRadians(float degrees) { return degrees * (float(3.14159265358979323846264338327950288) / float(180)); }

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor)
};
