"""SDM-PSI's parameter file, shared by the comparison scripts that invoke `sdm_parse`.

Only `nStudies` is dataset-specific; `nImputs` and `VoxelsMask` may be 0. Extracted here rather
than duplicated, because getting this file wrong makes `sdm_parse` fail in ways that look like a
data problem -- it cost hours to establish that `pp` does not write the config `mi` needs.
"""
PARAMS = """<MissSdm_parameters>
    <Global_parameters>
        <StudyMask>gray_matter_2mm</StudyMask>
        <CorrelationTemplate>gray_matter_2mm</CorrelationTemplate>
        <VoxelMM>2</VoxelMM>
        <nThreads>0</nThreads>
    </Global_parameters>
    <Model Name="MyMean">
        <Folder>analysis_MyMean</Folder>
        <nImputs>0</nImputs>
        <nPermutations>0</nPermutations>
        <PermsPath>permutations.asc</PermsPath>
        <VoxelMM>2</VoxelMM>
        <CorrelationTemplate>gray_matter_2mm</CorrelationTemplate>
        <VoxelsMask>0</VoxelsMask>
        <nThreads>1</nThreads>
        <nStudies>{n_studies}</nStudies>
        <useIntercept>1</useIntercept>
        <PermOnlyMeta>false</PermOnlyMeta>
        <SaveSubjectMaps>false</SaveSubjectMaps>
        <AccurateHeterogeneity>false</AccurateHeterogeneity>
        <MleCoefSmooth>1</MleCoefSmooth>
        <MleTau2Smooth>1</MleTau2Smooth>
        <MleTau2Scale>1</MleTau2Scale>
        <MleNLeaveOneOut>2</MleNLeaveOneOut>
        <VartMaxIterations>30</VartMaxIterations>
        <VartTol>0.0099999998</VartTol>
        <VartLearningScale>2</VartLearningScale>
        <Vars n="0"/>
        <Hypothesis n="5">
            <Hyp00000>1</Hyp00000>
            <Hyp00001>0</Hyp00001>
            <Hyp00002>0</Hyp00002>
            <Hyp00003>0</Hyp00003>
            <Hyp00004>0</Hyp00004>
        </Hypothesis>
        <Var1></Var1>
        <Var2></Var2>
        <Var3></Var3>
        <Var4></Var4>
        <Filter></Filter>
        <Hypothesis1>1</Hypothesis1>
        <Hypothesis2>0</Hypothesis2>
        <Hypothesis3>0</Hypothesis3>
        <Hypothesis4>0</Hypothesis4>
        <Hypothesis5>0</Hypothesis5>
        <Maps/>
        <Masks/>
    </Model>
</MissSdm_parameters>
"""
