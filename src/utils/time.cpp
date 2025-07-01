#include "time.hpp"
#include <ctime>
#include <string>

std::tm parseTime(const std::string& datetime) {
    std::tm tm = {};
    strptime(datetime.c_str(), "%Y-%m-%d %H:%M:%S", &tm);
    tm.tm_isdst = -1;
    return tm;
}

std::string formatTime(const std::tm& tm) {
    char buffer[30];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &tm);
    return std::string(buffer);
}

std::string addTimeDelta(const std::string& start, const std::string& calendar, int deltaMinutes) {
    std::tm tm = parseTime(start);

    if (calendar == "julian") {
        std::time_t t = std::mktime(&tm);
        t += deltaMinutes * 60;
        tm = *std::localtime(&t);
    } else if (calendar == "no_leap") {
        // Simple no-leap workaround (single-year)
        int minutesPerDay = 1440;
        int totalMinutes = tm.tm_yday * minutesPerDay + tm.tm_hour * 60 + tm.tm_min + deltaMinutes;

        static const int monthDays[12] = { 31,28,31,30,31,30,31,31,30,31,30,31 };

        int newDayOfYear = totalMinutes / minutesPerDay % 365;
        int newMinuteOfDay = totalMinutes % minutesPerDay;

        tm.tm_hour = newMinuteOfDay / 60;
        tm.tm_min = newMinuteOfDay % 60;
        tm.tm_sec = 0;

        int month = 0;
        int day = newDayOfYear;
        while (day >= monthDays[month]) {
            day -= monthDays[month];
            month++;
        }
        tm.tm_mon = month;
        tm.tm_mday = day + 1;
    } else {
        return "Invalid calendar type!";
    }

    return formatTime(tm);
}