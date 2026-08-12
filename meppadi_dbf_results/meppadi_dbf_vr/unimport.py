import unreal

# **********************************************************************

# The following parameters are partly set in r.avaflow, but should be revised

nfiles = 30 # number of time steps
blendfile = '' # prefix from blender file
interval = 4 # length of time step in frames
framerate = 60 # frame rate
addactor = '' # additional actor to appear at end ('' to skip)
addtime = 100 # time for additional actor to appear
addx = 0.0 # x position of additional actor
addy = 0.0 # y position of additional actor
addstartz = 0.0 # start z position of additional actor
addendz = 0.0 # end z position of additional actor
addscale = 1.0 # scale of additional actor

# **********************************************************************

# The following code should only be modified by experienced users

frate = unreal.FrameRate(numerator = framerate, denominator = 1)

lseq = unreal.LevelSequenceEditorBlueprintLibrary.get_current_level_sequence()
lseq.set_display_rate(frate)

lseq.set_playback_start(0)
lseq.set_playback_end(interval * nfiles + addtime)

actor_system = unreal.get_editor_subsystem(unreal.EditorActorSubsystem)
actors = unreal.EditorLevelLibrary.get_all_level_actors()

for i in range(0, nfiles+1):

    for actor in actors:
    
        if i<10: itext = '000' + str(i)
        elif i<100: itext = '00' + str(i)
        elif i<1000: itext = '0' + str(i)
        else: itext = str(i)
        
        if actor.get_actor_label() == blendfile + '_Surface' + itext:

            poss = lseq.add_possessable(actor)
            unreal.LevelSequenceEditorBlueprintLibrary.refresh_current_level_sequence()

            vistrack = poss.add_track(unreal.MovieSceneVisibilityTrack)
            vistrack.set_property_name_and_path('bHiddenInGame', 'bHiddenInGame')

            vissection = vistrack.add_section()
            vissection.set_start_frame_seconds(0)
            vissection.set_end_frame_seconds(interval * nfiles + addtime)

            vischannel = vissection.get_channels()[0]
            vischannel.add_key(time=unreal.FrameNumber(0), new_value=False)
            vischannel.add_key(time=unreal.FrameNumber(interval * i), new_value=True)
            if i < nfiles: vischannel.add_key(time=unreal.FrameNumber(interval * i + interval), new_value=False)

if not addactor == '':

    for actor in actors:

        if actor.get_actor_label() == addactor:

            poss = lseq.add_possessable(actor)
            unreal.LevelSequenceEditorBlueprintLibrary.refresh_current_level_sequence()

            vistrack = poss.add_track(unreal.MovieSceneVisibilityTrack)
            vistrack.set_property_name_and_path('bHiddenInGame', 'bHiddenInGame')

            vissection = vistrack.add_section()
            vissection.set_start_frame_seconds(0)
            vissection.set_end_frame_seconds(interval * nfiles + addtime)

            vischannel = vissection.get_channels()[0]
            vischannel.add_key(time=unreal.FrameNumber(interval * nfiles -1), new_value=False)
            vischannel.add_key(time=unreal.FrameNumber(interval * nfiles), new_value=True)

            transtrack = poss.add_track(unreal.MovieScene3DTransformTrack)

            transsection = transtrack.add_section()
            transsection.set_start_frame_seconds(0)
            transsection.set_end_frame_seconds(interval * nfiles + addtime)

            transchannel = transsection.get_channels()[0]
            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addx)

            transchannel = transsection.get_channels()[1]
            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addy)

            transchannel = transsection.get_channels()[2]
            transchannel.add_key(time=unreal.FrameNumber(interval * nfiles), new_value=addstartz)
            transchannel.add_key(time=unreal.FrameNumber(int(interval * nfiles + addtime * 0.9)), new_value=addendz)

            transchannel = transsection.get_channels()[6]
            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)
            
            transchannel = transsection.get_channels()[7]
            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)

            transchannel = transsection.get_channels()[8]
            transchannel.add_key(time=unreal.FrameNumber(0), new_value=addscale)
