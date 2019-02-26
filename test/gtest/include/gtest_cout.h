/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++.                                     |
|                                                                         |
|   Copyright(C) 2018  Alberto Cuoci                                      |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef _GTEST_COUT_H_

#define _GTEST_COUT_H_

#include "gtest/gtest.h"
#include <stdarg.h>

namespace testing
{
	namespace internal
	{

		#ifdef _WIN32

		enum GTestColor {
			COLOR_DEFAULT,
			COLOR_RED,
			COLOR_GREEN,
			COLOR_YELLOW
		};
		
		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
//		void ColoredPrintf(GTestColor color, const char* fmt, ...)
//		{};

		#elif _WIN64

		enum GTestColor {
			COLOR_DEFAULT,
			COLOR_RED,
			COLOR_GREEN,
			COLOR_YELLOW
		};

		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);

		#elif __apple__

		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);

		#elif __linux__

		// Returns the ANSI color code for the given color.  COLOR_DEFAULT is
		// an invalid input.
		static const char* MyGetAnsiColorCode(GTestColor color)
		{
			switch (color)
			{
				case COLOR_RED:     return "1";
				case COLOR_GREEN:   return "2";
				case COLOR_YELLOW:  return "3";
				default:            return NULL;
			};
		}

		// Returns true iff Google Test should use colors in the output.
		bool MyShouldUseColor(bool stdout_is_tty)
		{

			const char* const gtest_color = GTEST_FLAG(color).c_str();

			if (String::CaseInsensitiveCStringEquals(gtest_color, "auto"))
			{
			#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MINGW
				// On Windows the TERM variable is usually not set, but the
				// console there does support colors.
				return stdout_is_tty;
			#else
				// On non-Windows platforms, we rely on the TERM variable.
				const char* const term = posix::GetEnv("TERM");
				const bool term_supports_color =
				String::CStringEquals(term, "xterm") ||
				String::CStringEquals(term, "xterm-color") ||
				String::CStringEquals(term, "xterm-256color") ||
				String::CStringEquals(term, "screen") ||
				String::CStringEquals(term, "screen-256color") ||
				String::CStringEquals(term, "tmux") ||
				String::CStringEquals(term, "tmux-256color") ||
				String::CStringEquals(term, "rxvt-unicode") ||
				String::CStringEquals(term, "rxvt-unicode-256color") ||
				String::CStringEquals(term, "linux") ||
				String::CStringEquals(term, "cygwin");
				return stdout_is_tty && term_supports_color;
			#endif  // GTEST_OS_WINDOWS
			}

			return String::CaseInsensitiveCStringEquals(gtest_color, "yes") ||
			String::CaseInsensitiveCStringEquals(gtest_color, "true") ||
			String::CaseInsensitiveCStringEquals(gtest_color, "t") ||
			String::CStringEquals(gtest_color, "1");
			// We take "yes", "true", "t", and "1" as meaning "yes".  If the
			// value is neither one of these nor "auto", we treat it as "no" to
			// be conservative.
		}

		// Helpers for printing colored strings to stdout. Note that on Windows, we
		// cannot simply emit special characters and have the terminal change colors.
		// This routine must actually emit the characters rather than return a string
		// that would be colored when printed, as can be done on Linux.
		static void MyColoredPrintf(GTestColor color, const char* fmt, ...)
		{
			va_list args;
			va_start(args, fmt);

			#if GTEST_OS_WINDOWS_MOBILE || GTEST_OS_SYMBIAN || GTEST_OS_ZOS || \
			GTEST_OS_IOS || GTEST_OS_WINDOWS_PHONE || GTEST_OS_WINDOWS_RT
				const bool use_color = AlwaysFalse();
			#else
				static const bool in_color_mode =
				MyShouldUseColor(posix::IsATTY(posix::FileNo(stdout)) != 0);
				const bool use_color = in_color_mode && (color != COLOR_DEFAULT);
			#endif  // GTEST_OS_WINDOWS_MOBILE || GTEST_OS_SYMBIAN || GTEST_OS_ZOS
			// The '!= 0' comparison is necessary to satisfy MSVC 7.1.

			if (!use_color)
			{
				vprintf(fmt, args);
				va_end(args);
				return;
			}

			#if GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE && \
			!GTEST_OS_WINDOWS_PHONE && !GTEST_OS_WINDOWS_RT && !GTEST_OS_WINDOWS_MINGW
				const HANDLE stdout_handle = GetStdHandle(STD_OUTPUT_HANDLE);

				// Gets the current text color.
				CONSOLE_SCREEN_BUFFER_INFO buffer_info;
				GetConsoleScreenBufferInfo(stdout_handle, &buffer_info);
				const WORD old_color_attrs = buffer_info.wAttributes;
				const WORD new_color = GetNewColor(color, old_color_attrs);

				// We need to flush the stream buffers into the console before each
				// SetConsoleTextAttribute call lest it affect the text that is already
				// printed but has not yet reached the console.
				fflush(stdout);
				SetConsoleTextAttribute(stdout_handle, new_color);

				vprintf(fmt, args);

				fflush(stdout);
				// Restores the text color.
				SetConsoleTextAttribute(stdout_handle, old_color_attrs);
			#else
				printf("\033[0;3%sm", MyGetAnsiColorCode(color));
				vprintf(fmt, args);
				printf("\033[m");  // Resets the terminal to default.
			#endif  // GTEST_OS_WINDOWS && !GTEST_OS_WINDOWS_MOBILE

			va_end(args);
		}

		#endif
	}
}

#ifdef _WIN32

#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[ INFO     ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

#elif __apple__

#define PRINTF(...)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[ INFO     ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

#elif __linux__

#define PRINTF(...)  do { testing::internal::MyColoredPrintf(testing::internal::COLOR_GREEN, "[ INFO     ] "); testing::internal::MyColoredPrintf(testing::internal::COLOR_YELLOW, __VA_ARGS__); } while(0)

#endif

// C++ stream interface
class TestCout : public std::stringstream
{

	public:

		~TestCout()
		{
			PRINTF("%s", str().c_str());
		}

};

#define TEST_COUT  TestCout()

#endif
