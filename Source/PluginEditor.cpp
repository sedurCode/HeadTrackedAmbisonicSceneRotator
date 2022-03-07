/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor (HeadTrackedAmbisonicSceneRotatorAudioProcessor& p)
    : AudioProcessorEditor (&p), audioProcessor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
	setSize(600, 200);

	gainSlider.setSliderStyle(juce::Slider::LinearBarVertical);
	gainSlider.setRange(0.f, 1.f, 0.0001f); // TODO: Make log
	gainSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, true, 20, 20);
	gainSlider.setPopupDisplayEnabled(true, true, this);
	gainSlider.setTextValueSuffix("Gain");
	gainSlider.setDoubleClickReturnValue(true, 0.f);
	gainSlider.setValue(1.f);

	yawSlider.setSliderStyle(juce::Slider::LinearBarVertical);
	yawSlider.setRange(-180.f, 180.f, 0.1f); // TODO: Make log
	yawSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, true, 20, 20);
	yawSlider.setPopupDisplayEnabled(true, true, this);
	yawSlider.setTextValueSuffix("Yaw");
	yawSlider.setDoubleClickReturnValue(true, 0.f);
	yawSlider.setValue(0.f);

	pitchSlider.setSliderStyle(juce::Slider::LinearBarVertical);
	pitchSlider.setRange(-180.f, 180.f, 0.1f); // TODO: Make log
	pitchSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, true, 20, 20);
	pitchSlider.setPopupDisplayEnabled(true, true, this);
	pitchSlider.setTextValueSuffix("Pitch");
	pitchSlider.setDoubleClickReturnValue(true, 0.f);
	pitchSlider.setValue(0.f);

	rollSlider.setSliderStyle(juce::Slider::LinearBarVertical);
	rollSlider.setRange(-180.f, 180.f, 0.1f); // TODO: Make log
	rollSlider.setTextBoxStyle(juce::Slider::TextBoxBelow, true, 20, 20);
	rollSlider.setPopupDisplayEnabled(true, true, this);
	rollSlider.setTextValueSuffix("Roll");
	rollSlider.setDoubleClickReturnValue(true, 0.f);
	rollSlider.setValue(0.f);

	addAndMakeVisible(gainSlider);
	addAndMakeVisible(yawSlider);
	addAndMakeVisible(pitchSlider);
	addAndMakeVisible(rollSlider);

	gainSlider.addListener(this);
	yawSlider.addListener(this);
	pitchSlider.addListener(this);
	rollSlider.addListener(this);

	headPanel.setListener(this);
	headPanel.setTopLeftPosition(8, 8);
	addAndMakeVisible(headPanel);

	V.x = 0.25;
	V.y = 0.25;
	V.z = 0.25;
	//setSize(headPanel.getWidth() + 16, headPanel.getHeight() + 16);
}

HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::~HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor()
{
}

//==============================================================================
void HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    g.setColour (juce::Colours::black);
    g.setFont (15.0f);
    g.drawFittedText ("Head tracked ambisonic panner", 100, 10,400, 20, juce::Justification::centred, 1);
}

void HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::resized()
{
    // This is generally where you'll want to lay out the positions of any
    // subcomponents in your editor..
	gainSlider.setBounds(200, 30, 40, getHeight() - 60);
	yawSlider.setBounds(300, 30, 40, getHeight() - 60);
	pitchSlider.setBounds(400, 30, 40, getHeight() - 60);
	rollSlider.setBounds(500, 30, 40, getHeight() - 60);
}

void HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::sliderValueChanged(juce::Slider* slider)
{
	if (slider == &gainSlider)
	{
		audioProcessor.setGain(slider->getValue());
	}
	if (slider == &yawSlider)
	{
		audioProcessor.setYaw(slider->getValue());
	}
	if (slider == &pitchSlider)
	{
		audioProcessor.setPitch(slider->getValue());
	}
	if (slider == &rollSlider)
	{
		audioProcessor.setRoll(slider->getValue());
	}

}

void HeadTrackedAmbisonicSceneRotatorAudioProcessorEditor::trackerChanged(const HeadMatrix& headMatrix)
{
	// headMatrix.transform and headMatrix.transformTranspose can be used here
	// to rotate an object.
	const juce::MessageManagerLock mmessageLock;
	float y, p, r;
	headMatrix.getEulerAngles(y, p, r);			// the reference 
	  yawSlider.setValue(y);
	pitchSlider.setValue(p);
	 rollSlider.setValue(r);

}