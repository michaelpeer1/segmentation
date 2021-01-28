#!/usr/bin/env python
import os


def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return template, outtype, annotation_classes


# Create keys
t1w = create_key(
    '{subject}/{session}/anat/{subject}_{session}_T1w')
func_localizer_run1 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-func_localizer_run-01_bold')
func_localizer_run2 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-func_localizer_run-02_bold')
object_viewing_run1_d1 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-object_viewing_run-01_day-01_bold')
object_viewing_run2_d1 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-object_viewing_run-02_day-01_bold')
object_viewing_run1_d2 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-object_viewing_run-01_day-02_bold')
object_viewing_run2_d2 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-object_viewing_run-02_day-02_bold')
jrd_run1 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-jrd_run-01_bold')
jrd_run2 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-jrd_run-02_bold')
jrd_run3 = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-jrd_run-03_bold')
resting_state = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-resting_state_bold')
navigation = create_key(
    '{subject}/{session}/func/{subject}_{session}_task-navigation_bold')
b0_phase = create_key(
    '{subject}/{session}/fmap/{subject}_{session}_phasediff')
b0_mag = create_key(
    '{subject}/{session}/fmap/{subject}_{session}_magnitude{item}')

def infotodict(seqinfo):

    last_run = len(seqinfo)

    info = {t1w:[],
        func_localizer_run1:[],
        func_localizer_run2:[],
        object_viewing_run1_d1:[],
        object_viewing_run2_d1:[],
        object_viewing_run1_d2:[],
        object_viewing_run2_d2:[],
        jrd_run1:[],
        jrd_run2:[],
        jrd_run3:[],
        resting_state:[],
        navigation:[],
        b0_phase:[], b0_mag:[]}

    for s in seqinfo:
        protocol = s.protocol_name.lower()

        if "mprage" in protocol:
            info[t1w].append(s.series_id)
        elif "localizer_1" in protocol:
            info[func_localizer_run1].append(s.series_id)
        elif "localizer_2" in protocol:
            info[func_localizer_run2].append(s.series_id)
        elif "objectviewing_1_d1" in protocol:
            info[object_viewing_run1_d1].append(s.series_id)
        elif "objectviewing_2_d1" in protocol:
            info[object_viewing_run2_d1].append(s.series_id)
        elif "objectviewing_1_d2" in protocol:
            info[object_viewing_run1_d2].append(s.series_id)
        elif "objectviewing_2_d2" in protocol:
            info[object_viewing_run2_d2].append(s.series_id)
        elif "jrd_1" in protocol:
            info[jrd_run1].append(s.series_id)
        elif "jrd_2" in protocol:
            info[jrd_run2].append(s.series_id)
        elif "jrd_3" in protocol:
            info[jrd_run3].append(s.series_id)
        elif "resting" in protocol:
            info[resting_state].append(s.series_id)
        elif "navigation" in protocol:
            info[navigation].append(s.series_id)
        elif "b0map" in protocol and "P" in s.image_type:
            info[b0_phase].append(s.series_id)
        elif "b0map" in protocol and "M" in s.image_type:
            info[b0_mag].append(s.series_id)
 
    return info

MetadataExtras = {
    b0_phase: {
        "EchoTime1": 0.004,
        "EchoTime2": 0.006
    }
}

IntendedFor = {
    b0_phase: [
    '{session}/func/sub-{subject}_{session}_task-func_localizer_run-01_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-func_localizer_run-02_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-object_viewing_run-01_day-01_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-object_viewing_run-02_day-01_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-object_viewing_run-01_day-02_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-object_viewing_run-02_day-02_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-jrd_run-01_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-jrd_run-02_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-jrd_run-03_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-resting_state_bold.nii.gz',
    '{session}/func/sub-{subject}_{session}_task-navigation_bold.nii.gz']
}

def ReplaceSubject(str_input):
    import re
    out = re.sub("sub-", "", str_input)
    return(out)

