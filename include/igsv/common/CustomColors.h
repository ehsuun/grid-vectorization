//=============================================================================
//
//  STRUCT : CustomColors
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== TYPEDEFS ===============================================================

  using Color3d = Eigen::RowVector3d;

  //== STRUCT DEFINITION ======================================================

  struct CustomColors {
    static const Color3d U();
    static const Color3d V();
    // static const Color3d Random();
    static const Color3d Palette(const int k) {
      switch (k % 25) {
      case 0: return Red;               // rgb( 255,   0,   0)
      case 1: return Orange;            // rgb( 255, 165,   0)
      case 2: return Gold;              // rgb( 255, 215,   0)
      case 3: return Green;             // rgb(   0, 128,   0)
      case 4: return Blue;              // rgb(   0,   0, 255)
      case 5: return Purple;            // rgb( 128,   0, 128)
      case 6: return Pink;              // rgb( 255, 192, 203)
      case 7: return DarkSalmon;        // rgb( 233, 150, 122)
      case 8: return Indigo;            // rgb(  75,   0, 130)
      case 9: return MediumVioletRed;   // rgb( 199,  21, 133)
      case 10: return DarkKhaki;        // rgb( 189, 183, 107)
      case 11: return RosyBrown;        // rgb( 188, 143, 143)
      case 12: return MediumAquamarine; // rgb( 102, 205, 170)
      case 13: return LightSkyBlue;     // rgb( 135, 206, 235)
      case 14: return Cyan;             // rgb(   0, 255, 255)
      case 15: return SandyBrown;       // rgb( 244, 164,  96)
      case 16: return BurlyWood;        // rgb( 222, 184, 135)
      case 17: return SlateBlue;        // rgb(  72,  61, 139)
      case 18: return SlateGray;        // rgb( 112, 128, 144)
      case 19: return DarkKhaki;        // rgb( 189, 183, 107)
      case 20: return SeaGreen;         // rgb(  46, 139,  87)
      case 21: return DarkRed;          // rgb( 139,   0,   0)
      case 22: return DeepPink;         // rgb( 255,  20, 147)
      case 23: return Aquamarine;       // rgb( 127, 255, 212)
      case 24: return DeepSkyBlue;      // rgb(   0, 191, 255)
      default: return Black;            // this should not be reached.
      }
    }
    //
    // https://en.wikipedia.org/wiki/Web_colors#X11_color_names
    //
    // Pink colors
    static const Color3d Pink;
    static const Color3d LightPink;
    static const Color3d HotPink;
    static const Color3d DeepPink;
    static const Color3d PaleVioletRed;
    static const Color3d MediumVioletRed;
    // Red colors
    static const Color3d LightSalmon;
    static const Color3d Salmon;
    static const Color3d DarkSalmon;
    static const Color3d LightCoral;
    static const Color3d IndianRed;
    static const Color3d Crimson;
    static const Color3d FireBrick;
    static const Color3d DarkRed;
    static const Color3d Red;
    // Orange colors
    static const Color3d OrangeRed;
    static const Color3d Tomato;
    static const Color3d Coral;
    static const Color3d DarkOrange;
    static const Color3d Orange;
    // Yellow colors
    static const Color3d Yellow;
    static const Color3d LightYellow;
    static const Color3d LemonChiffon;
    static const Color3d LightGoldenrodYellow;
    static const Color3d PapayaWhip;
    static const Color3d Moccasin;
    static const Color3d PeachPuff;
    static const Color3d PaleGoldenrod;
    static const Color3d Khaki;
    static const Color3d DarkKhaki;
    static const Color3d Gold;
    // Brown colors
    static const Color3d Cornsilk;
    static const Color3d BlanchedAlmond;
    static const Color3d Bisque;
    static const Color3d NavajoWhite;
    static const Color3d Wheat;
    static const Color3d BurlyWood;
    static const Color3d Tan;
    static const Color3d RosyBrown;
    static const Color3d SandyBrown;
    static const Color3d Goldenrod;
    static const Color3d DarkGoldenrod;
    static const Color3d Peru;
    static const Color3d Chocolate;
    static const Color3d SaddleBrown;
    static const Color3d Sienna;
    static const Color3d Brown;
    static const Color3d Maroon;
    // Green colors
    static const Color3d DarkOliveGreen;
    static const Color3d Olive;
    static const Color3d OliveDrab;
    static const Color3d YellowGreen;
    static const Color3d LimeGreen;
    static const Color3d Lime;
    static const Color3d LawnGreen;
    static const Color3d Chartreuse;
    static const Color3d GreenYellow;
    static const Color3d SpringGreen;
    static const Color3d MediumSpringGreen;
    static const Color3d LightGreen;
    static const Color3d PaleGreen;
    static const Color3d DarkSeaGreen;
    static const Color3d MediumAquamarine;
    static const Color3d MediumSeaGreen;
    static const Color3d SeaGreen;
    static const Color3d ForestGreen;
    static const Color3d Green;
    static const Color3d DarkGreen;
    // Cyan colors
    static const Color3d Aqua;
    static const Color3d Cyan;
    static const Color3d LightCyan;
    static const Color3d PaleTurquoise;
    static const Color3d Aquamarine;
    static const Color3d Turquoise;
    static const Color3d MediumTurquoise;
    static const Color3d DarkTurquoise;
    static const Color3d LightSeaGreen;
    static const Color3d CadetBlue;
    static const Color3d DarkCyan;
    static const Color3d Teal;
    // Blue colors
    static const Color3d LightSteelBlue;
    static const Color3d PowderBlue;
    static const Color3d LightBlue;
    static const Color3d SkyBlue;
    static const Color3d LightSkyBlue;
    static const Color3d DeepSkyBlue;
    static const Color3d DodgerBlue;
    static const Color3d CornflowerBlue;
    static const Color3d SteelBlue;
    static const Color3d RoyalBlue;
    static const Color3d Blue;
    static const Color3d MediumBlue;
    static const Color3d DarkBlue;
    static const Color3d Navy;
    static const Color3d MidnightBlue;
    // Purple, violet, and magenta colors
    static const Color3d Lavender;
    static const Color3d Thistle;
    static const Color3d Plum;
    static const Color3d Violet;
    static const Color3d Orchid;
    static const Color3d Fuchsia;
    static const Color3d Magenta;
    static const Color3d MediumOrchid;
    static const Color3d MediumPurple;
    static const Color3d BlueViolet;
    static const Color3d DarkViolet;
    static const Color3d DarkOrchid;
    static const Color3d DarkMagenta;
    static const Color3d Purple;
    static const Color3d Indigo;
    static const Color3d DarkSlateBlue;
    static const Color3d SlateBlue;
    static const Color3d MediumSlateBlue;
    // White colors
    static const Color3d White;
    static const Color3d Snow;
    static const Color3d Honeydew;
    static const Color3d MintCream;
    static const Color3d Azure;
    static const Color3d AliceBlue;
    static const Color3d GhostWhite;
    static const Color3d WhiteSmoke;
    static const Color3d Seashell;
    static const Color3d Beige;
    static const Color3d OldLace;
    static const Color3d FloralWhite;
    static const Color3d Ivory;
    static const Color3d AntiqueWhite;
    static const Color3d Linen;
    static const Color3d LavenderBlush;
    static const Color3d MistyRose;
    // Gray and black colors
    static const Color3d Gainsboro;
    static const Color3d LightGray;
    static const Color3d Silver;
    static const Color3d DarkGray;
    static const Color3d Gray;
    static const Color3d DimGray;
    static const Color3d LightSlateGray;
    static const Color3d SlateGray;
    static const Color3d DarkSlateGray;
    static const Color3d Black;
  };

} // namespace IGSV
