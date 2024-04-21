#include "c74_min.h"
#include <algorithm>
#include <vector>
#include <unordered_map>

using namespace c74::min;
using namespace std;

class VoiceAllocator : public object<VoiceAllocator> {
private:
    // Define the Note struct
    struct Note {
        long ID;
        long velocity;
        long midiNote;
        long voicesAssigned;
        long voicesRequired;
        vector<long> voiceNumbers;
    };

    // Define member variables
    vector<Note> notes;
    unordered_map<long, Note*> voiceToNoteMap;

    long calculateVoicesPerNote() {
        // Calculate the total number of active notes
        long totalActiveNotes = notes.size();

        // Calculate the minimum number of voices per note
        long minVoicesPerNote = 8 / totalActiveNotes;

        // Calculate the number of remaining voices
        long remainingVoices = 8 % totalActiveNotes;

        // Determine the number of voices for each note
        for (size_t i = 0; i < notes.size(); ++i) {
            notes[i].voicesRequired = minVoicesPerNote;
        }

        // Distribute remaining voices to newer notes
        for (size_t i = 0; i < remainingVoices; ++i) {
            notes[notes.size() - i - 1].voicesRequired = minVoicesPerNote + 1;
        }
        

        // Return the total number of voices required for the new note
        return notes[notes.size() - 1].voicesRequired;
    }

    void allocateVoices(Note& note, long numVoices) {
        // Assign voices to the note
        long assignedVoices = note.voicesAssigned; // Get the current number of assigned voices
        long voiceNumber = 1; // Start from the first voice

        while (assignedVoices < numVoices && voiceNumber <= 8) {
            // Check if the voice is already assigned
            if (voiceToNoteMap.find(voiceNumber) == voiceToNoteMap.end()) {
                // Voice is free, assign it to the note
                note.voiceNumbers.push_back(voiceNumber);
                voiceToNoteMap[voiceNumber] = &note;

                // Construct output message and send it
                m_outlet.send(note.ID, note.midiNote, note.velocity, voiceNumber);

                assignedVoices++; // Increment the count of assigned voices
            }

            voiceNumber++; // Move to the next voice
        }

        // Update the count of voices assigned to the note
        note.voicesAssigned = assignedVoices;
    }



    void freeVoices(Note& note) {
        // Update the state of the note
        for (long voiceNumber : note.voiceNumbers) {
            voiceToNoteMap.erase(voiceNumber);

            // Construct output message and send it
            m_outlet.send(note.ID, note.midiNote, 0, voiceNumber);
        }

        // Clear the note from the notes vector
        auto it = find_if(notes.begin(), notes.end(), [&note](const Note& n) { return n.ID == note.ID; });
        if (it != notes.end()) {
            notes.erase(it);
        }
    }


    void stealVoices(Note& newNote, long totalVoicesRequired) {
        // Sort notes by ID (oldest to newest)
        sort(notes.begin(), notes.end(), [](const Note& a, const Note& b) {
            return a.ID < b.ID;
        });

        // Keep track of how many voices are still needed to be stolen
        long remainingVoicesToSteal = totalVoicesRequired;

        // Iterate through each note
        for (auto& note : notes) {
            // Check if the note has assigned voices
            if (note.voicesAssigned > 0) {
                // Determine how many voices to steal from the current note
                long voicesAvailable = note.voicesAssigned - note.voicesRequired;
                if (voicesAvailable < 0)
                    voicesAvailable = 0;
                long voicesToSteal = min(voicesAvailable, remainingVoicesToSteal);

                // Iterate through the voices to steal
                for (int i = 0; i < voicesToSteal; ++i) {
                    // Get the last voice number assigned to the note
                    long stolenVoiceNumber = note.voiceNumbers.back();

                    // Remove the stolen voice from the note's voice numbers
                    note.voiceNumbers.pop_back();

                    // Remove the stolen voice from the voice-to-note map
                    voiceToNoteMap.erase(stolenVoiceNumber);

                    // Send a note-off message for the stolen voice
                    m_outlet.send(note.ID, note.midiNote, 0, stolenVoiceNumber);
                }

                // Update the remaining voices to steal
                remainingVoicesToSteal -= voicesToSteal;
                
                // Update the voices assigned to the current note
                note.voicesAssigned -= voicesToSteal;

                // Check if enough voices have been stolen
                if (remainingVoicesToSteal <= 0) break;
            }
        }

        // Allocate stolen voices to the new note
        allocateVoices(newNote, totalVoicesRequired - remainingVoicesToSteal);
    }




public:
    // Constructor, inlets, outlets, and process method...

    MIN_DESCRIPTION    { "8-Voice unison allocation" };
    MIN_TAGS        { "polyphony" };
    MIN_AUTHOR        { "Osprey Instruments, Peter Corbett" };
    MIN_RELATED        { "poly, poly~" };
    
    inlet<> m_inlet { this, "(anything) Input" };
    outlet<> m_outlet { this, "(anything) Output" };    // message outlet
    outlet<> d_outlet { this, "(anything) Output" };    // debug outlet
    
    // Process the input list
    message<> list { this, "list",
        MIN_FUNCTION {
            // Check if the input list has the correct number of elements
            if (args.size() == 3)
            {
                // Extract input elements
                long ID = args[0];
                long midiNote = args[1];
                long velocity = args[2];
                
                // Find the corresponding note or create a new one
                Note* currentNote;
                auto it = find_if(notes.begin(), notes.end(), [ID](const Note& n) { return n.ID == ID; });
                if (it != notes.end()) {
                    currentNote = &(*it);
                } else {
                    notes.push_back({ID, velocity, midiNote, 0, 0, {}});
                    currentNote = &notes.back();
                }
                
                // Handle note off event
                if (velocity == 0) {
                    freeVoices(*currentNote);
                    
                    // Check if all notes are released
                    bool allNotesReleased = all_of(notes.begin(), notes.end(), [](const Note& n) { return n.voicesAssigned == 0; });
                    
                    // If all notes are released, reset the voice allocation
                    if (allNotesReleased) {
                        for (auto& note : notes) {
                            note.voiceNumbers.clear();
                        }
                        voiceToNoteMap.clear();
                    }
                } else { // Handle note on event
                    long totalVoicesForNewNote = calculateVoicesPerNote();
                    if (totalVoicesForNewNote > 0) {
                        // Allocate voices to the new note
                        allocateVoices(*currentNote, totalVoicesForNewNote);
                        
                        // If not enough voices available, steal from older notes
                        if (currentNote->voicesAssigned < totalVoicesForNewNote) {
                            long remainingVoices = totalVoicesForNewNote - currentNote->voicesAssigned;
                            stealVoices(*currentNote, remainingVoices);
                        }
                    }
                }
                
                
                // debug routines
                for (size_t i = 0; i < notes.size(); ++i) {
                    
                    // Send the required and assigned voice allocation
                    d_outlet.send("Note: ", i + 1, "Voices Required: ", notes[i].voicesRequired, "Allocated: ", notes[i].voicesAssigned);
                    
                    // Send the note ID and voice numbers
                    for (size_t j = 0; j < notes[i].voiceNumbers.size(); ++j) {
                        d_outlet.send("Note: ", notes[i].ID, "Voice Number: ", notes[i].voiceNumbers[j]);
                    }
                }
            }
            else
            {
                d_outlet.send("Incoming list has incorrect number of elements.");
            }
            return {};
        }
    };
};

MIN_EXTERNAL(VoiceAllocator);
