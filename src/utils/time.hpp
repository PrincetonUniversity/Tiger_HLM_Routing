#pragma once

#include <string>
#include <ctime>

// Parses a datetime string "YYYY-MM-DD HH:MM:SS" into a std::tm struct
std::tm parseTime(const std::string& datetime);

// Formats a std::tm as "YYYY-MM-DD HH:MM:SS"
std::string formatTime(const std::tm& tm);

// Adds a delta in minutes to the start time using "julian" or "no_leap" calendar
std::string addTimeDelta(const std::string& start, const std::string& calendar, int deltaMinutes);