// Abf2Tmp.cpp : Written by Ryan Thompson 2008
//  You may copy, distribute and modify this code in any way you wish
//  please visit the project page for questions and comments
//  http://sourceforge.net/projects/matlababffilere

#include "stdafx.h"
#include <windows.h>
#include "Abffiles.h"
#include <iostream>
#include <stdio.h>

#define BUFSIZE 16384

bool writeMetaData(const char* fname,ABFFileHeader& FH);
bool writeShorts(FILE * pf, const char* name, short* data, int total);
bool writeFloats(FILE * pf, const char* name, float* data, int total);
bool writeChars(FILE * pf, const char* name, char* data, int total);
bool writeStrings(FILE * pf, const char* name, char* data, int total, int maxStringSize);
bool writeLongs(FILE * pf, const char* name, long* data, int total);


BOOL ShowABFError( char *szFileName, int nError )
{
char szTxt[80];
if (!ABF_BuildErrorText( nError, szFileName, szTxt, sizeof(szTxt) ))
sprintf( szTxt, "Unknown error number: %d\r\n", nError );
printf( "ERROR: %s\n", szTxt );
return FALSE;
}

bool UnicodeToAnsi(const _TCHAR* s, char* out, int len =256)
{
	int i =0;
	while(s[i] !=0 && i<len)
	{
		out[i] = CHAR(s[i]);
		i++;	
	}
	out[i] =0;
	if (i >=len) return false;
	return true;
}


int _tmain(int argc, _TCHAR* argv[])
{
	if (argc <2)
	{
		std::cout<< "Usage: abf2Tmp inFile.abf\n";
		exit(0);
	}
	char inFname[1024];
	if(!UnicodeToAnsi(argv[1], inFname, 1024))
	{
		std::cout << "filename too long.\n";
		exit(0);
	}
	

	int hFile;
	int nError;
	ABFFileHeader FH;
	DWORD dwMaxEpi = 0;
	UINT uMaxSamples = BUFSIZE;
	char outFileName[1024];
	sprintf(outFileName,"%s.txt",inFname);
	//LPCSTR szFileName = CHAR(outFileName);
	BOOL bSuccess= ABF_ReadOpen(inFname, &hFile, ABF_DATAFILE, &FH,
		&uMaxSamples, &dwMaxEpi, &nError);
	if (!bSuccess) return ShowABFError(inFname, nError);

	//float* pfBuffer = new float[FH.nADCNumChannels*BUFSIZE];
	float* pfBuffer = new float[FH.nADCNumChannels*uMaxSamples];
	
	DWORD dwEpisode;
	FILE* pf = fopen(outFileName,"w");
	UINT puNumSamples = uMaxSamples;
	for (dwEpisode=1;dwEpisode <= dwMaxEpi; dwEpisode++)
	{
		for (int nChannel =0; nChannel < FH.nADCNumChannels; nChannel++)
		{
			//some channels might not be used this compensates for that.
			int thisChannel = FH.nADCSamplingSeq[nChannel];
			bSuccess = ABF_ReadChannel( hFile, &FH, thisChannel, dwEpisode, &pfBuffer[nChannel*uMaxSamples],
				&puNumSamples, &nError );
			if (!bSuccess) 
				{
/*					if (nError == 1012 ){ //channel not found
						//here we should fill the channel with 0's
						for (int iData = nChannel*BUFSIZE;  iData < (nChannel+1) * BUFSIZE; iData ++)
						{
							pfBuffer[iData] = 0;
						}
					}
					else*/
						return ShowABFError(inFname, nError);
				}
		}
		for (unsigned int i=0; i < puNumSamples ;i++)
		{
			for (int nChannel =0; nChannel < FH.nADCNumChannels; nChannel ++)
			{
				fprintf(pf,"%f%c", pfBuffer[nChannel * uMaxSamples + i], nChannel != FH.nADCNumChannels -1?'\t':'\n');
			}
		}


	}
	fclose(pf);

	sprintf(outFileName,"%s-md.txt",inFname);
	if (!writeMetaData(outFileName,FH))
	{
		std::cout << "Couldn't write metadata file\n";
	}

	if (!ABF_Close( hFile, &nError ))
		return ShowABFError(inFname, nError);

	return 0;
}


bool writeMetaData(const char* fname, ABFFileHeader& FH)
{
	FILE* pf = fopen(fname,"w");
	if (!pf) return false;

	// GROUP #1 - File ID and size information
	fprintf(pf,"fFileVersionNumber\t%f\n",FH.fFileVersionNumber);
	fprintf(pf,"nOperationMode\t%d\n",FH.nOperationMode);
	fprintf(pf,"lActualAcqLength\t%d\n",FH.lActualAcqLength);
	fprintf(pf,"nNumPointsIgnored\t%d\n",FH.nNumPointsIgnored);
	fprintf(pf,"lActualEpisodes\t%d\n",FH.lActualEpisodes);
	fprintf(pf,"uFileStartDate\t%d\n",FH.uFileStartDate);
	fprintf(pf,"uFileStartTimeMS\t%d\n",FH.uFileStartTimeMS);
	fprintf(pf,"lStopwatchTime\t%d\n",FH.lStopwatchTime);
	fprintf(pf,"fHeaderVersionNumber\t%f\n",FH.fHeaderVersionNumber);
	fprintf(pf,"nFileType\t%d\n",FH.nFileType);


   // GROUP #2 - File Structure
	fprintf(pf,"lDataSectionPtr\t%d\n",FH.lDataSectionPtr);
	fprintf(pf,"lTagSectionPtr\t%d\n",FH.lTagSectionPtr);
	fprintf(pf,"lNumTagEntries\t%d\n",FH.lNumTagEntries);
	fprintf(pf,"lScopeConfigPtr\t%d\n",FH.lScopeConfigPtr);
	fprintf(pf,"lNumScopes\t%d\n",FH.lNumScopes);
	fprintf(pf,"lDeltaArrayPtr\t%d\n",FH.lDeltaArrayPtr);
	fprintf(pf,"lNumDeltas\t%d\n",FH.lNumDeltas);
	fprintf(pf,"lVoiceTagPtr\t%d\n",FH.lVoiceTagPtr);
	fprintf(pf,"lVoiceTagEntries\t%d\n",FH.lVoiceTagEntries);
	fprintf(pf,"lSynchArrayPtr\t%d\n",FH.lSynchArrayPtr);
	fprintf(pf,"lSynchArraySize\t%d\n",FH.lSynchArraySize);
	fprintf(pf,"nDataFormat\t%d\n",FH.nDataFormat);
	fprintf(pf,"nSimultaneousScan\t%d\n",FH.nSimultaneousScan);
	fprintf(pf,"lStatisticsConfigPtr\t%d\n",FH.lStatisticsConfigPtr);
	fprintf(pf,"lAnnotationSectionPtr\t%d\n",FH.lAnnotationSectionPtr);
	fprintf(pf,"lNumAnnotations\t%d\n",FH.lNumAnnotations);

/*
TODO : I don't think this metadata is important
   long     lDACFilePtr[ABF_DACCOUNT];
   long     lDACFileNumEpisodes[ABF_DACCOUNT];
*/


   // GROUP #3 - Trial hierarchy information
	fprintf(pf,"nADCNumChannels\t%d\n",FH.nADCNumChannels);   
	fprintf(pf,"fADCSequenceInterval\t%f\n",FH.fADCSequenceInterval);   
	fprintf(pf,"uFileCompressionRatio\t%d\n",FH.uFileCompressionRatio);   
	fprintf(pf,"bEnableFileCompression\t%d\n",FH.bEnableFileCompression);   
	fprintf(pf,"fSynchTimeUnit\t%f\n",FH.fSynchTimeUnit);   
	fprintf(pf,"fSecondsPerRun\t%f\n",FH.fSecondsPerRun);   
	fprintf(pf,"lNumSamplesPerEpisode\t%d\n",FH.lNumSamplesPerEpisode);   
	fprintf(pf,"lPreTriggerSamples\t%d\n",FH.lPreTriggerSamples);   
	fprintf(pf,"lEpisodesPerRun\t%d\n",FH.lEpisodesPerRun);   
	fprintf(pf,"lRunsPerTrial\t%d\n",FH.lRunsPerTrial);   
	fprintf(pf,"lNumberOfTrials\t%d\n",FH.lNumberOfTrials);   
	fprintf(pf,"nAveragingMode\t%d\n",FH.nAveragingMode);   
	fprintf(pf,"nUndoRunCount\t%d\n",FH.nUndoRunCount);   
	fprintf(pf,"nFirstEpisodeInRun\t%d\n",FH.nFirstEpisodeInRun);   
	fprintf(pf,"fTriggerThreshold\t%f\n",FH.fTriggerThreshold);   
	fprintf(pf,"nTriggerSource\t%d\n",FH.nTriggerSource);   
	fprintf(pf,"nTriggerAction\t%d\n",FH.nTriggerAction);   
	fprintf(pf,"nTriggerPolarity\t%d\n",FH.nTriggerPolarity);   
	fprintf(pf,"fScopeOutputInterval\t%f\n",FH.fScopeOutputInterval);   
	fprintf(pf,"fEpisodeStartToStart\t%f\n",FH.fEpisodeStartToStart);   
	fprintf(pf,"fRunStartToStart\t%f\n",FH.fRunStartToStart);   
	fprintf(pf,"fTrialStartToStart\t%f\n",FH.fTrialStartToStart);   
	fprintf(pf,"lAverageCount\t%d\n",FH.lAverageCount);   
	fprintf(pf,"nAutoTriggerStrategy\t%d\n",FH.nAutoTriggerStrategy);   
	fprintf(pf,"fFirstRunDelayS\t%f\n",FH.fFirstRunDelayS);   
   
   

   // GROUP #4 - Display Parameters
	fprintf(pf,"nDataDisplayMode\t%d\n",FH.nDataDisplayMode);   
	fprintf(pf,"nChannelStatsStrategy\t%d\n",FH.nChannelStatsStrategy);   
	fprintf(pf,"lSamplesPerTrace\t%d\n",FH.lSamplesPerTrace);   
	fprintf(pf,"lStartDisplayNum\t%d\n",FH.lStartDisplayNum);   
	fprintf(pf,"lFinishDisplayNum\t%d\n",FH.lFinishDisplayNum);   
	fprintf(pf,"nShowPNRawData\t%d\n",FH.nShowPNRawData);   
	fprintf(pf,"fStatisticsPeriod\t%f\n",FH.fStatisticsPeriod);   
	fprintf(pf,"lStatisticsMeasurements\t%d\n",FH.lStatisticsMeasurements);   
	fprintf(pf,"nStatisticsSaveStrategy\t%d\n",FH.nStatisticsSaveStrategy);   

   // GROUP #5 - Hardware information
	fprintf(pf,"fADCRange\t%f\n",FH.fADCRange);   
	fprintf(pf,"fDACRange\t%f\n",FH.fDACRange);   
	fprintf(pf,"lADCResolution\t%d\n",FH.lADCResolution);   
	fprintf(pf,"lDACResolution\t%d\n",FH.lDACResolution);   
	fprintf(pf,"nDigitizerADCs\t%d\n",FH.nDigitizerADCs);   
	fprintf(pf,"nDigitizerDACs\t%d\n",FH.nDigitizerDACs);   
	fprintf(pf,"nDigitizerTotalDigitalOuts\t%d\n",FH.nDigitizerTotalDigitalOuts);   
	fprintf(pf,"nDigitizerSynchDigitalOuts\t%d\n",FH.nDigitizerSynchDigitalOuts);   
	fprintf(pf,"nDigitizerType\t%d\n",FH.nDigitizerType);   

   // GROUP #6 Environmental Information
	fprintf(pf,"nExperimentType\t%d\n",FH.nExperimentType);   
	fprintf(pf,"nManualInfoStrategy\t%d\n",FH.nManualInfoStrategy);   
	fprintf(pf,"fCellID1\t%d\n",FH.fCellID1);   
	fprintf(pf,"fCellID2\t%d\n",FH.fCellID2);   
	fprintf(pf,"fCellID3\t%d\n",FH.fCellID3);   
	fprintf(pf,"sProtocolPath\t%s\n",FH.sProtocolPath);   
	fprintf(pf,"sCreatorInfo\t%s\n",FH.sCreatorInfo);   
	fprintf(pf,"sModifierInfo\t%s\n",FH.sModifierInfo);   
	fprintf(pf,"nCommentsEnable\t%d\n",FH.nCommentsEnable);   
	fprintf(pf,"sFileComment\t%s\n",FH.sFileComment);   
/*
	TODO: I am not sure what type of data this is, if anyone knows send me email.
   short    nTelegraphEnable[ABF_ADCCOUNT];
   short    nTelegraphInstrument[ABF_ADCCOUNT];
   float    fTelegraphAdditGain[ABF_ADCCOUNT];
   float    fTelegraphFilter[ABF_ADCCOUNT];
   float    fTelegraphMembraneCap[ABF_ADCCOUNT];
   float    fTelegraphAccessResistance[ABF_ADCCOUNT];
   short    nTelegraphMode[ABF_ADCCOUNT];
   short    nTelegraphDACScaleFactorEnable[ABF_DACCOUNT];
	
*/
/*
TODO : I don't believe this data is relevant.
   short    nAutoAnalyseEnable;

   GUID     FileGUID;
   float    fInstrumentHoldingLevel[ABF_DACCOUNT];
   unsigned long ulFileCRC;
   short    nCRCEnable;
*/

   // GROUP #7 - Multi-channel information
	fprintf(pf,"nSignalType\t%d\n",FH.nSignalType);   
	writeShorts(pf,"nADCPtoLChannelMap",FH.nADCPtoLChannelMap,ABF_ADCCOUNT);
	writeShorts(pf,"nADCSamplingSeq",FH.nADCSamplingSeq,ABF_ADCCOUNT);
	writeFloats(pf,"fADCProgrammableGain",FH.fADCProgrammableGain,ABF_ADCCOUNT);
	writeFloats(pf,"fADCDisplayAmplification",FH.fADCDisplayAmplification,ABF_ADCCOUNT);
	writeFloats(pf,"fADCDisplayOffset",FH.fADCDisplayOffset,ABF_ADCCOUNT);
	writeFloats(pf,"fInstrumentScaleFactor",FH.fInstrumentScaleFactor,ABF_ADCCOUNT);
	writeFloats(pf,"fInstrumentOffset",FH.fInstrumentOffset,ABF_ADCCOUNT);
	writeFloats(pf,"fSignalGain",FH.fSignalGain,ABF_ADCCOUNT);
	writeFloats(pf,"fSignalOffset",FH.fSignalOffset,ABF_ADCCOUNT);
	writeFloats(pf,"fSignalLowpassFilter",FH.fSignalLowpassFilter,ABF_ADCCOUNT);
	writeFloats(pf,"fSignalHighpassFilter",FH.fSignalHighpassFilter,ABF_ADCCOUNT);
	writeChars(pf,"nLowpassFilterType",FH.nLowpassFilterType,ABF_ADCCOUNT);
	writeChars(pf,"nHighpassFilterType",FH.nHighpassFilterType,ABF_ADCCOUNT);
	
	writeStrings(pf,"sADCChannelName",(char*)FH.sADCChannelName,ABF_ADCCOUNT,ABF_ADCNAMELEN);
	writeStrings(pf,"sADCUnits",(char*)FH.sADCUnits,ABF_ADCCOUNT,ABF_ADCNAMELEN);
	writeFloats(pf,"fDACScaleFactor",FH.fDACScaleFactor,ABF_ADCCOUNT);
	writeFloats(pf,"fDACHoldingLevel",FH.fDACHoldingLevel,ABF_ADCCOUNT);
	writeFloats(pf,"fDACCalibrationFactor",FH.fDACCalibrationFactor,ABF_ADCCOUNT);
	writeFloats(pf,"fDACCalibrationOffset",FH.fDACCalibrationOffset,ABF_ADCCOUNT);
	writeStrings(pf,"sDACChannelName",(char*)FH.sDACChannelName,ABF_DACCOUNT,ABF_DACNAMELEN);
	writeStrings(pf,"sDACChannelUnits",(char*)FH.sDACChannelUnits,ABF_DACCOUNT,ABF_DACNAMELEN);
 /*
 TODO : not sure what this group is
   // GROUP #9 - Epoch Waveform and Pulses
   short    nDigitalEnable;
   short    nActiveDACChannel;                     // should retire !
   short    nDigitalDACChannel;
   short    nDigitalHolding;
   short    nDigitalInterEpisode;
   short    nDigitalTrainActiveLogic;                                   
   short    nDigitalValue[ABF_EPOCHCOUNT];
   short    nDigitalTrainValue[ABF_EPOCHCOUNT];                         
   bool     bEpochCompression[ABF_EPOCHCOUNT];
   short    nWaveformEnable[ABF_DACCOUNT];
   short    nWaveformSource[ABF_DACCOUNT];
   short    nInterEpisodeLevel[ABF_DACCOUNT];
   short    nEpochType[ABF_DACCOUNT][ABF_EPOCHCOUNT];
   float    fEpochInitLevel[ABF_DACCOUNT][ABF_EPOCHCOUNT];
   float    fEpochLevelInc[ABF_DACCOUNT][ABF_EPOCHCOUNT];
   long     lEpochInitDuration[ABF_DACCOUNT][ABF_EPOCHCOUNT];
   long     lEpochDurationInc[ABF_DACCOUNT][ABF_EPOCHCOUNT];
*/
   // GROUP #10 - DAC Output File
	writeFloats(pf,"fDACFileScale",FH.fDACFileScale,ABF_DACCOUNT);
	writeFloats(pf,"fDACFileOffset",FH.fDACFileOffset,ABF_DACCOUNT);
	writeLongs(pf,"lDACFileEpisodeNum",FH.lDACFileEpisodeNum,ABF_DACCOUNT);
	writeShorts(pf,"nDACFileADCNum",FH.nDACFileADCNum,ABF_DACCOUNT);
	writeStrings(pf,"sDACFilePath",(char*)FH.sDACFilePath,ABF_DACCOUNT,ABF_PATHLEN);


   // GROUP #11 - Presweep (conditioning) pulse train
	writeShorts(pf,"nConditEnable",FH.nConditEnable,ABF_DACCOUNT);
	writeLongs(pf,"lConditNumPulses",FH.lConditNumPulses,ABF_DACCOUNT);
	writeFloats(pf,"fBaselineDuration",FH.fBaselineDuration,ABF_DACCOUNT);
	writeFloats(pf,"fBaselineLevel",FH.fBaselineLevel,ABF_DACCOUNT);
	writeFloats(pf,"fStepDuration",FH.fStepDuration,ABF_DACCOUNT);
	writeFloats(pf,"fStepLevel",FH.fStepLevel,ABF_DACCOUNT);
	writeFloats(pf,"fPostTrainPeriod",FH.fPostTrainPeriod,ABF_DACCOUNT);
	writeFloats(pf,"fPostTrainLevel",FH.fPostTrainLevel,ABF_DACCOUNT);
	writeShorts(pf,"nMembTestEnable",FH.nMembTestEnable,ABF_DACCOUNT);
	writeFloats(pf,"fMembTestPreSettlingTimeMS",FH.fMembTestPreSettlingTimeMS,ABF_DACCOUNT);
	writeFloats(pf,"fMembTestPostSettlingTimeMS",FH.fMembTestPostSettlingTimeMS,ABF_DACCOUNT);

   // GROUP #12 - Variable parameter user list
	writeShorts(pf,"nULEnable",FH.nULEnable,ABF_USERLISTCOUNT);
	writeShorts(pf,"nULParamToVary",FH.nULParamToVary,ABF_USERLISTCOUNT);
	writeShorts(pf,"nULRepeat",FH.nULRepeat,ABF_USERLISTCOUNT);
	writeStrings(pf,"sDACFilePath",(char*)FH.sDACFilePath,ABF_USERLISTCOUNT,ABF_USERLISTLEN);

   // GROUP #13 - Statistics measurements
	fprintf(pf,"nStatsEnable\t%d\n",FH.nStatsEnable);   
	fprintf(pf,"nStatsActiveChannels\t%d\n",FH.nStatsActiveChannels);   
	fprintf(pf,"nStatsSearchRegionFlags\t%d\n",FH.nStatsSearchRegionFlags);   
	fprintf(pf,"nStatsSmoothing\t%d\n",FH.nStatsSmoothing);   
	fprintf(pf,"nStatsSmoothingEnable\t%d\n",FH.nStatsSmoothingEnable);   
	fprintf(pf,"nStatsBaseline\t%d\n",FH.nStatsBaseline);   
	fprintf(pf,"nStatsBaselineDAC\t%d\n",FH.nStatsBaselineDAC);   
	fprintf(pf,"lStatsBaselineStart\t%d\n",FH.lStatsBaselineStart);   
	fprintf(pf,"lStatsBaselineEnd\t%d\n",FH.lStatsBaselineEnd);   
	writeLongs(pf,"lStatsMeasurements",FH.lStatsMeasurements,ABF_STATS_REGIONS);
	writeLongs(pf,"lStatsStart",FH.lStatsStart,ABF_STATS_REGIONS);
	writeLongs(pf,"lStatsEnd",FH.lStatsEnd,ABF_STATS_REGIONS);
	writeShorts(pf,"nRiseBottomPercentile",FH.nRiseBottomPercentile,ABF_STATS_REGIONS);
	writeShorts(pf,"nRiseTopPercentile",FH.nRiseTopPercentile,ABF_STATS_REGIONS);
	writeShorts(pf,"nDecayBottomPercentile",FH.nDecayBottomPercentile,ABF_STATS_REGIONS);
	writeShorts(pf,"nDecayTopPercentile",FH.nDecayTopPercentile,ABF_STATS_REGIONS);
	writeShorts(pf,"nStatsChannelPolarity",FH.nStatsChannelPolarity,ABF_ADCCOUNT);
	writeShorts(pf,"nStatsSearchMode",FH.nStatsSearchMode,ABF_STATS_REGIONS);
	writeShorts(pf,"nStatsSearchDAC",FH.nStatsSearchDAC,ABF_STATS_REGIONS);

   // GROUP #14 - Channel Arithmetic
	fprintf(pf,"nArithmeticEnable\t%d\n",FH.nArithmeticEnable);   
	fprintf(pf,"nArithmeticExpression\t%d\n",FH.nArithmeticExpression);   
	fprintf(pf,"fArithmeticLowerLimit\t%f\n",FH.fArithmeticLowerLimit);   
	fprintf(pf,"fArithmeticLowerLimit\t%f\n",FH.fArithmeticLowerLimit);   
	fprintf(pf,"nArithmeticADCNumA\t%d\n",FH.nArithmeticADCNumA);   
	fprintf(pf,"nArithmeticADCNumB\t%d\n",FH.nArithmeticADCNumB);   
	fprintf(pf,"fArithmeticK1\t%f\n",FH.fArithmeticK1);   
	fprintf(pf,"fArithmeticK2\t%f\n",FH.fArithmeticK2);   
	fprintf(pf,"fArithmeticK3\t%f\n",FH.fArithmeticK3);   
	fprintf(pf,"fArithmeticK4\t%f\n",FH.fArithmeticK4);   
	fprintf(pf,"fArithmeticK5\t%f\n",FH.fArithmeticK5);   
	fprintf(pf,"sArithmeticOperator\t%s\n",FH.sArithmeticOperator);   //Test for overflow condition?
	fprintf(pf,"sArithmeticUnits\t%s\n",FH.sArithmeticUnits);   
   // GROUP #15 - Leak subtraction
	fprintf(pf,"nPNPosition\t%d\n",FH.nPNPosition);   
	fprintf(pf,"nPNNumPulses\t%d\n",FH.nPNNumPulses);   
	fprintf(pf,"nPNPolarity\t%d\n",FH.nPNPolarity);   
	fprintf(pf,"fPNSettlingTime\t%f\n",FH.fPNSettlingTime);   
	fprintf(pf,"fPNInterpulse\t%f\n",FH.fPNInterpulse);   
	writeShorts(pf,"nLeakSubtractType",FH.nLeakSubtractType,ABF_DACCOUNT);
	writeFloats(pf,"fPNHoldingLevel",FH.fPNHoldingLevel,ABF_DACCOUNT);
	writeShorts(pf,"nLeakSubtractADCIndex",FH.nLeakSubtractADCIndex,ABF_DACCOUNT);

   // GROUP #16 - Miscellaneous variables
	fprintf(pf,"nLevelHysteresis\t%d\n",FH.nLevelHysteresis);   
	fprintf(pf,"lTimeHysteresis\t%d\n",FH.lTimeHysteresis);   
	fprintf(pf,"nAllowExternalTags\t%d\n",FH.nAllowExternalTags);   
	fprintf(pf,"nAverageAlgorithm\t%d\n",FH.nAverageAlgorithm);   
	fprintf(pf,"fAverageWeighting\t%f\n",FH.fAverageWeighting);   
	fprintf(pf,"nUndoPromptStrategy\t%d\n",FH.nUndoPromptStrategy);   
	fprintf(pf,"nTrialTriggerSource\t%d\n",FH.nTrialTriggerSource);   
	fprintf(pf,"nStatisticsDisplayStrategy\t%d\n",FH.nStatisticsDisplayStrategy);   
	fprintf(pf,"nExternalTagType\t%d\n",FH.nExternalTagType);   
	fprintf(pf,"lHeaderSize\t%d\n",FH.lHeaderSize);   
	fprintf(pf,"nStatisticsClearStrategy\t%d\n",FH.nStatisticsClearStrategy);   
   
   // GROUP #17 - Trains parameters
/*
TODO : find out about this
	long     lEpochPulsePeriod[ABF_DACCOUNT][ABF_EPOCHCOUNT];
   long     lEpochPulseWidth [ABF_DACCOUNT][ABF_EPOCHCOUNT];
*/
   // GROUP #18 - Application version data
	fprintf(pf,"nCreatorMajorVersion\t%d\n",FH.nCreatorMajorVersion);   
	fprintf(pf,"nCreatorMinorVersion\t%d\n",FH.nCreatorMinorVersion);   
	fprintf(pf,"nCreatorBugfixVersion\t%d\n",FH.nCreatorBugfixVersion);   
	fprintf(pf,"nCreatorBuildVersion\t%d\n",FH.nCreatorBuildVersion);   
	fprintf(pf,"nModifierMajorVersion\t%d\n",FH.nModifierMajorVersion);   
	fprintf(pf,"nModifierMinorVersion\t%d\n",FH.nModifierMinorVersion);   
	fprintf(pf,"nModifierBugfixVersion\t%d\n",FH.nModifierBugfixVersion);   
	fprintf(pf,"nModifierBuildVersion\t%d\n",FH.nModifierBuildVersion);   

   // GROUP #19 - LTP protocol
	fprintf(pf,"nLTPType\t%d\n",FH.nLTPType);   
	writeShorts(pf,"nLTPUsageOfDAC",FH.nLTPUsageOfDAC,ABF_DACCOUNT);
	writeShorts(pf,"nLTPPresynapticPulses",FH.nLTPPresynapticPulses,ABF_DACCOUNT);

   // GROUP #20 - Digidata 132x Trigger out flag
	fprintf(pf,"nScopeTriggerOut\t%d\n",FH.nScopeTriggerOut);   


   // GROUP #22 - Alternating episodic mode
	fprintf(pf,"nAlternateDACOutputState\t%d\n",FH.nAlternateDACOutputState);   
	fprintf(pf,"nAlternateDigitalOutputState\t%d\n",FH.nAlternateDigitalOutputState);   
	writeShorts(pf,"nAlternateDigitalValue",FH.nAlternateDigitalValue,ABF_EPOCHCOUNT);
	writeShorts(pf,"nAlternateDigitalTrainValue",FH.nAlternateDigitalTrainValue,ABF_EPOCHCOUNT);

   // GROUP #23 - Post-processing actions
	writeFloats(pf,"fPostProcessLowpassFilter",FH.fPostProcessLowpassFilter,ABF_ADCCOUNT);
	writeChars(pf,"nPostProcessLowpassFilterType",FH.nPostProcessLowpassFilterType,ABF_ADCCOUNT);

   // GROUP #24 - Legacy gear shift info
	fprintf(pf,"fLegacyADCSequenceInterval\t%f\n",FH.fLegacyADCSequenceInterval);   
	fprintf(pf,"fLegacyADCSecondSequenceInterval\t%f\n",FH.fLegacyADCSecondSequenceInterval);   
	fprintf(pf,"lLegacyClockChange\t%d\n",FH.lLegacyClockChange);   
	fprintf(pf,"lLegacyNumSamplesPerEpisode\t%d\n",FH.lLegacyNumSamplesPerEpisode);   

   
   
   fclose(pf);
	return true;
}

bool writeLongs(FILE * pf, const char* name, long* data, int total)
{
	fprintf(pf,"%s\t", name);
	for (int i=0; i< total; i++)
	{
		fprintf(pf,"%d%c", data[i], i==total-1?'\n':'\t');
	}
	return true;
}


bool writeShorts(FILE * pf, const char* name, short* data, int total)
{
	fprintf(pf,"%s\t", name);
	for (int i=0; i< total; i++)
	{
		fprintf(pf,"%d%c", data[i], i==total-1?'\n':'\t');
	}
	return true;
}

bool writeFloats(FILE * pf, const char* name, float* data, int total)
{
	fprintf(pf,"%s\t", name);
	for (int i=0; i< total; i++)
	{
		fprintf(pf,"%f%c", data[i], i==total-1?'\n':'\t');
	}
	return true;
}

bool writeChars(FILE * pf, const char* name, char* data, int total)
{
	fprintf(pf,"%s\t", name);
	for (int i=0; i< total; i++)
	{
		fprintf(pf,"%d%c", data[i], i==total-1?'\n':'\t');
	}
	return true;
}

bool writeStrings(FILE * pf, const char* name, char* data, int total, int maxStringSize)
{
	fprintf(pf,"%s\t", name);
	for (int i=0; i< total*maxStringSize; i+= maxStringSize)
	{
		char buf[512];
		for (int j = 0; j < maxStringSize; j ++)
		{
			buf[j] = data[i+j];
		}
		buf[maxStringSize] =0;
		fprintf(pf,"%s%c", buf, i==(total-1)*maxStringSize?'\n':'\t');
	}
	return true;

}