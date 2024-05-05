/*!
 * d3-hexgrid plugin v0.3.1. https://github.com/ligoo/d3-hexgrid.git.
 */

(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (global = typeof globalThis !== 'undefined' ? globalThis : global || self, factory(global.d3 = {}));
})(this, (function (exports) { 'use strict';

  function ascending(a, b) {
    return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
  }

  function bisector(compare) {
    if (compare.length === 1) compare = ascendingComparator(compare);
    return {
      left: function(a, x, lo, hi) {
        if (lo == null) lo = 0;
        if (hi == null) hi = a.length;
        while (lo < hi) {
          var mid = lo + hi >>> 1;
          if (compare(a[mid], x) < 0) lo = mid + 1;
          else hi = mid;
        }
        return lo;
      },
      right: function(a, x, lo, hi) {
        if (lo == null) lo = 0;
        if (hi == null) hi = a.length;
        while (lo < hi) {
          var mid = lo + hi >>> 1;
          if (compare(a[mid], x) > 0) hi = mid;
          else lo = mid + 1;
        }
        return lo;
      }
    };
  }

  function ascendingComparator(f) {
    return function(d, x) {
      return ascending(f(d), x);
    };
  }

  bisector(ascending);

  function number(x) {
    return x === null ? NaN : +x;
  }

  function extent(values, valueof) {
    var n = values.length,
        i = -1,
        value,
        min,
        max;

    if (valueof == null) {
      while (++i < n) { // Find the first comparable value.
        if ((value = values[i]) != null && value >= value) {
          min = max = value;
          while (++i < n) { // Compare the remaining values.
            if ((value = values[i]) != null) {
              if (min > value) min = value;
              if (max < value) max = value;
            }
          }
        }
      }
    }

    else {
      while (++i < n) { // Find the first comparable value.
        if ((value = valueof(values[i], i, values)) != null && value >= value) {
          min = max = value;
          while (++i < n) { // Compare the remaining values.
            if ((value = valueof(values[i], i, values)) != null) {
              if (min > value) min = value;
              if (max < value) max = value;
            }
          }
        }
      }
    }

    return [min, max];
  }

  function quantile(values, p, valueof) {
    if (valueof == null) valueof = number;
    if (!(n = values.length)) return;
    if ((p = +p) <= 0 || n < 2) return +valueof(values[0], 0, values);
    if (p >= 1) return +valueof(values[n - 1], n - 1, values);
    var n,
        i = (n - 1) * p,
        i0 = Math.floor(i),
        value0 = +valueof(values[i0], i0, values),
        value1 = +valueof(values[i0 + 1], i0 + 1, values);
    return value0 + (value1 - value0) * (i - i0);
  }

  /* eslint-disable no-param-reassign */

  /**
   * Checks and if required converts the 1D extent to a 2D extent.
   * @param {Array} userExtent  Either the full 2D extent or just width and height.
   * @return                    The full 2D extent.
   */
  function expandExtent(userExtent) {
    const nestedArrayLength = Array.from(new Set(userExtent.map(e => e.length)))[0];
    let extentLong = Array(2);
    if (nestedArrayLength === 2) {
      extentLong = userExtent;
    } else if (nestedArrayLength === undefined) {
      extentLong = [[0, 0], userExtent];
    } else {
      throw new Error("Check 'extent' is in the anticipated form [[x0,y0],[x1,y1]] or [x1,y1]");
    }
    return extentLong;
  }

  /**
   * Checks and sets given value to greater than 0.
   * @param  {number} v       Value.
   * @param  {string} name    Value name.
   * @return {number}         Value.
   */
  function convertToMin(v, name, min) {
    if (v >= min) {
      return v;
    }
    // eslint-disable-next-line no-console
    console.warn(`${name} should be ${min} or greater. Coerced to ${min}.`);
    return min;
  }

  /**
   * Produce corner points for a pointy hexagon.
   * @param  {Object} center Hexagon center position.
   * @param  {number} r      Radius of hexagon.
   * @param  {number} i      Index of point to calculate.
   * @return {Object}        Hexagon corner position.
   */
  function pointyHexCorner(center, r, i) {
    const point = {};
    const angleDegree = 60 * i - 30;
    const angleRadian = Math.PI / 180 * angleDegree;
    point.x = center.x + r * Math.cos(angleRadian);
    point.y = center.y + r * Math.sin(angleRadian);
    return point;
  }

  /**
   * Draw a hexagon.
   * @param  {Object} context The canvas context.
   * @param  {Object} corners Hexagon corner positions.
   * @param  {String} action  'fill' or 'stroke'.
   * @param  {String} colour  Colour.
   * @return {[type]}         undefined
   */
  function hexDraw(context, corners, colour) {
    let action = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 'fill';
    context.beginPath();
    corners.forEach(d => {
      d === 0 ? context.moveTo(d.x, d.y) : context.lineTo(d.x, d.y);
    });
    context.closePath();
    if (action === 'fill') {
      context.fillStyle = colour;
      context.fill();
    } else if (action === 'stroke') {
      context.strokeStyle = colour;
      context.stroke();
    } else {
      throw new Error("hexDraw action needs to be either 'fill' or 'stroke'");
    }
  }

  /**
   * Calculates the circle radius in pixel, given a circle polygon.
   * @param   {Object}    geoCirclePolygon  The circle polygon.
   * @param   {function}  projection        The D3 projection function.
   * @return  {number}                      The radius in pixel.
   */
  function getPixelRadius(geoCirclePolygon, projection) {
    // Get radius in pixel.
    const circleDataGeo = geoCirclePolygon.coordinates[0];
    const circleDataY = circleDataGeo.map(d => projection(d)[1]);
    const circleDiameter = extent(circleDataY);
    const radiusPixel = (circleDiameter[1] - circleDiameter[0]) / 2;
    return radiusPixel;
  }

  var epsilon = 1e-6;
  var pi = Math.PI;
  var halfPi = pi / 2;
  var tau = pi * 2;

  var degrees = 180 / pi;
  var radians = pi / 180;

  var abs = Math.abs;
  var atan2 = Math.atan2;
  var cos = Math.cos;
  var sin = Math.sin;
  var sqrt = Math.sqrt;

  function acos(x) {
    return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
  }

  function asin(x) {
    return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
  }

  function spherical(cartesian) {
    return [atan2(cartesian[1], cartesian[0]), asin(cartesian[2])];
  }

  function cartesian(spherical) {
    var lambda = spherical[0], phi = spherical[1], cosPhi = cos(phi);
    return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
  }

  // TODO return d
  function cartesianNormalizeInPlace(d) {
    var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    d[0] /= l, d[1] /= l, d[2] /= l;
  }

  function constant(x) {
    return function() {
      return x;
    };
  }

  function compose(a, b) {

    function compose(x, y) {
      return x = a(x, y), b(x[0], x[1]);
    }

    if (a.invert && b.invert) compose.invert = function(x, y) {
      return x = b.invert(x, y), x && a.invert(x[0], x[1]);
    };

    return compose;
  }

  function rotationIdentity(lambda, phi) {
    return [abs(lambda) > pi ? lambda + Math.round(-lambda / tau) * tau : lambda, phi];
  }

  rotationIdentity.invert = rotationIdentity;

  function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
    return (deltaLambda %= tau) ? (deltaPhi || deltaGamma ? compose(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma))
      : rotationLambda(deltaLambda))
      : (deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma)
      : rotationIdentity);
  }

  function forwardRotationLambda(deltaLambda) {
    return function(lambda, phi) {
      return lambda += deltaLambda, [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
    };
  }

  function rotationLambda(deltaLambda) {
    var rotation = forwardRotationLambda(deltaLambda);
    rotation.invert = forwardRotationLambda(-deltaLambda);
    return rotation;
  }

  function rotationPhiGamma(deltaPhi, deltaGamma) {
    var cosDeltaPhi = cos(deltaPhi),
        sinDeltaPhi = sin(deltaPhi),
        cosDeltaGamma = cos(deltaGamma),
        sinDeltaGamma = sin(deltaGamma);

    function rotation(lambda, phi) {
      var cosPhi = cos(phi),
          x = cos(lambda) * cosPhi,
          y = sin(lambda) * cosPhi,
          z = sin(phi),
          k = z * cosDeltaPhi + x * sinDeltaPhi;
      return [
        atan2(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
        asin(k * cosDeltaGamma + y * sinDeltaGamma)
      ];
    }

    rotation.invert = function(lambda, phi) {
      var cosPhi = cos(phi),
          x = cos(lambda) * cosPhi,
          y = sin(lambda) * cosPhi,
          z = sin(phi),
          k = z * cosDeltaGamma - y * sinDeltaGamma;
      return [
        atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi),
        asin(k * cosDeltaPhi - x * sinDeltaPhi)
      ];
    };

    return rotation;
  }

  // Generates a circle centered at [0°, 0°], with a given radius and precision.
  function circleStream(stream, radius, delta, direction, t0, t1) {
    if (!delta) return;
    var cosRadius = cos(radius),
        sinRadius = sin(radius),
        step = direction * delta;
    if (t0 == null) {
      t0 = radius + direction * tau;
      t1 = radius - step / 2;
    } else {
      t0 = circleRadius(cosRadius, t0);
      t1 = circleRadius(cosRadius, t1);
      if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * tau;
    }
    for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
      point = spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
      stream.point(point[0], point[1]);
    }
  }

  // Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].
  function circleRadius(cosRadius, point) {
    point = cartesian(point), point[0] -= cosRadius;
    cartesianNormalizeInPlace(point);
    var radius = acos(-point[1]);
    return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
  }

  function geoCircle() {
    var center = constant([0, 0]),
        radius = constant(90),
        precision = constant(6),
        ring,
        rotate,
        stream = {point: point};

    function point(x, y) {
      ring.push(x = rotate(x, y));
      x[0] *= degrees, x[1] *= degrees;
    }

    function circle() {
      var c = center.apply(this, arguments),
          r = radius.apply(this, arguments) * radians,
          p = precision.apply(this, arguments) * radians;
      ring = [];
      rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
      circleStream(stream, r, p, 1);
      c = {type: "Polygon", coordinates: [ring]};
      ring = rotate = null;
      return c;
    }

    circle.center = function(_) {
      return arguments.length ? (center = typeof _ === "function" ? _ : constant([+_[0], +_[1]]), circle) : center;
    };

    circle.radius = function(_) {
      return arguments.length ? (radius = typeof _ === "function" ? _ : constant(+_), circle) : radius;
    };

    circle.precision = function(_) {
      return arguments.length ? (precision = typeof _ === "function" ? _ : constant(+_), circle) : precision;
    };

    return circle;
  }

  /**
   * Sniffs the unit and converts to either "m" or "km".
   * @param   {string} unit The user given unit.
   * @return  {string}      The clean unit string.
   */
  function getUnitString(unit) {
    let unitLower = unit.toLowerCase();
    if (unitLower === 'm' || unitLower === 'km') {
      return unitLower;
    }
    if (unitLower === 'kilometres' || unitLower === 'kilometre' || unitLower === 'kilometers' || unitLower === 'kilometer') {
      unitLower = 'km';
    } else if (unitLower === 'miles' || unitLower === 'mile') {
      unitLower = 'm';
    } else {
      throw new Error('Please provide the unit identifier as either "km" for kilometres or "m" for miles');
    }
    return unitLower;
  }

  /**
   *
   * @param   {number}    radiusDistance  The user given distance in either miles or km.
   * @param   {string}    distanceUnit    The user chosen distance unit (miles or km).
   * @param   {function}  projection      The D3 projection function.
   * @param   {Array}     center          The center coordinates of the drawing area.
   * @return  {Object}                    The geo circle, the radius in degrees and in pixel.
   */
  function convertUnitRadius (radiusDistance, distanceUnit, projection, center) {
    // Get radius in degrees
    const unit = getUnitString(distanceUnit);
    const RADIUS_EARTH = unit === 'm' ? 3959 : 6371;
    const radiusRadians = radiusDistance / RADIUS_EARTH;
    const radiusDegrees = radiusRadians * (180 / Math.PI);

    // Get geo circle data.
    const circlePolygon = geoCircle().radius(radiusDegrees).center(projection.invert(center));

    // Get radius in pixel.
    const radiusPixel = getPixelRadius(circlePolygon(), projection);
    return {
      circlePolygon,
      radiusDegrees,
      radiusPixel
    };
  }

  var thirdPi = Math.PI / 3,
      angles = [0, thirdPi, 2 * thirdPi, 3 * thirdPi, 4 * thirdPi, 5 * thirdPi];

  function pointX(d) {
    return d[0];
  }

  function pointY(d) {
    return d[1];
  }

  function hexbin() {
    var x0 = 0,
        y0 = 0,
        x1 = 1,
        y1 = 1,
        x = pointX,
        y = pointY,
        r,
        dx,
        dy;

    function hexbin(points) {
      var binsById = {}, bins = [], i, n = points.length;

      for (i = 0; i < n; ++i) {
        if (isNaN(px = +x.call(null, point = points[i], i, points))
            || isNaN(py = +y.call(null, point, i, points))) continue;

        var point,
            px,
            py,
            pj = Math.round(py = py / dy),
            pi = Math.round(px = px / dx - (pj & 1) / 2),
            py1 = py - pj;

        if (Math.abs(py1) * 3 > 1) {
          var px1 = px - pi,
              pi2 = pi + (px < pi ? -1 : 1) / 2,
              pj2 = pj + (py < pj ? -1 : 1),
              px2 = px - pi2,
              py2 = py - pj2;
          if (px1 * px1 + py1 * py1 > px2 * px2 + py2 * py2) pi = pi2 + (pj & 1 ? 1 : -1) / 2, pj = pj2;
        }

        var id = pi + "-" + pj, bin = binsById[id];
        if (bin) bin.push(point);
        else {
          bins.push(bin = binsById[id] = [point]);
          bin.x = (pi + (pj & 1) / 2) * dx;
          bin.y = pj * dy;
        }
      }

      return bins;
    }

    function hexagon(radius) {
      var x0 = 0, y0 = 0;
      return angles.map(function(angle) {
        var x1 = Math.sin(angle) * radius,
            y1 = -Math.cos(angle) * radius,
            dx = x1 - x0,
            dy = y1 - y0;
        x0 = x1, y0 = y1;
        return [dx, dy];
      });
    }

    hexbin.hexagon = function(radius) {
      return "m" + hexagon(radius == null ? r : +radius).join("l") + "z";
    };

    hexbin.centers = function() {
      var centers = [],
          j = Math.round(y0 / dy),
          i = Math.round(x0 / dx);
      for (var y = j * dy; y < y1 + r; y += dy, ++j) {
        for (var x = i * dx + (j & 1) * dx / 2; x < x1 + dx / 2; x += dx) {
          centers.push([x, y]);
        }
      }
      return centers;
    };

    hexbin.mesh = function() {
      var fragment = hexagon(r).slice(0, 4).join("l");
      return hexbin.centers().map(function(p) { return "M" + p + "m" + fragment; }).join("");
    };

    hexbin.x = function(_) {
      return arguments.length ? (x = _, hexbin) : x;
    };

    hexbin.y = function(_) {
      return arguments.length ? (y = _, hexbin) : y;
    };

    hexbin.radius = function(_) {
      return arguments.length ? (r = +_, dx = r * 2 * Math.sin(thirdPi), dy = r * 1.5, hexbin) : r;
    };

    hexbin.size = function(_) {
      return arguments.length ? (x0 = y0 = 0, x1 = +_[0], y1 = +_[1], hexbin) : [x1 - x0, y1 - y0];
    };

    hexbin.extent = function(_) {
      return arguments.length ? (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1], hexbin) : [[x0, y0], [x1, y1]];
    };

    return hexbin.radius(1);
  }

  /**
   * Configure the hexbin generator.
   * @param  {Array}    extent   Drawing area extent.
   * @param  {number}   radius   The desired hex radius.
   * @return {function}          Hexbin generator function.
   */
  function setHexGenerator (extent, radius) {
    // Set the hexbin generator. Note, x and y will
    // be set later when prepping the user data.
    // Also round radius to the nearest 0.5 step.
    return hexbin().extent(extent).radius(radius).x(d => d.x).y(d => d.y);
  }

  /**
   * Checks which pixel of the image are filled
   * returning pixel positions to draw hexes on.
   * @param  {Array}    size        Width and height of base element.
   * @param  {function} pathGen     D3 path generator function.
   * @param  {Object}   geo         GeoJSON representing the object to tesselate.
   * @param  {number}   r           Hexagon radius.
   * @param  {string}   action      Drawing action `fill` or `stroke`.
   * @param  {number}   band        Extension of image (factor of r).
   * @return {Uint8ClampedArray}  Array of A values (from RGBA) per pixel.
   */
  // export default function(size, precision, pathGen, geo, r, action, band) {
  function getImageData (size, pathGen, geo, r, action, band) {
    const gridExtentStroke = band * r;
    const edgeBand = gridExtentStroke + 2 * r;

    // For debugging; append the canvas to the body and just draw on it.
    const canvas = document.createElement('canvas');
    // const canvas = d3.select('body').append('canvas').node();
    [canvas.width, canvas.height] = size;
    const context = canvas.getContext('2d', {
      willReadFrequently: true
    });
    const canvasPath = pathGen.context(context);

    // Draw.
    context.beginPath();
    canvasPath(geo);
    if (action === 'fill') {
      // debugger
      if (band) {
        context.lineWidth = gridExtentStroke;
        context.stroke();
      }
      context.fill();
    } else if (action === 'stroke') {
      context.lineWidth = edgeBand;
      context.stroke();
    }

    // Remove side effect of setting the path's context.
    pathGen.context(null);

    // Get the pixel rgba data but only keep the 4th value (alpha).
    const imgData = context.getImageData(0, 0, size[0], size[1]).data;
    return imgData.filter((d, i) => i % 4 === 3);
  }

  /**
   * Checks for each center if it covers a pixel in the image.
   * Checks only for centers that are within the bounds of width and height.
   * Note, this can be optimised (for loop instead of filter all).
   * @param  {Array}              centers     Hexagon centers covering the
   *                                          extent of the drawing canvas.
   * @param  {Uint8ClampedArray}  image       Pixel alpha values indicating fill.
   * @param  {Array}              size        Width and height of drawing canvas.
   * @param  {number}             precision   Hidden canvas ratio of the
   *                                          drawing canvas.
   * @return {Array}                          Hexagon centers covering
   *                                          the displayed object.
   */
  // export default function(centers, image, size, precision) {
  function getImageCenters (centers, image, size) {
    const [w, h] = size;
    return centers.filter(center => {
      return (
        // Guarantee centers to be within bounds.
        center[0] >= 0 && center[0] <= w && center[1] >= 0 && center[1] <= h && image[Math.floor(center[0]) + Math.floor(center[1]) * Math.floor(w)]
      );
    }).map((center, i) => {
      return {
        id: i,
        x: center[0],
        y: center[1],
        gridpoint: 1,
        cover: 1
      };
    });
  }

  /**
   * Checks for each center if it covers a pixel in the image.
   * @param  {Array}              centers     Hexagon centers covering the
   *                                          breadth of the drawing canvas.
   * @param  {Uint8ClampedArray}  image       Pixels indicating fill.
   * @param  {Array}              size        Width and height of drawing canvas.
   * @param  {number}             precision   Hidden canvas ratio of the
   *                                          drawing canvas.
   * @return {Array}                          Hexagon centers covering the
   *                                          displayed object only.
   */
  // export default function(centers, image, size, precision) {
  function getEdgeCenters (centers, image, size) {
    const w = size[0];
    return centers.filter(el => image[Math.floor(el.x) + w * Math.floor(el.y)]);
  }

  /**
   * Produe the canvas image of a specifically sized and scaled hexagon,
   * the canvas image of the desired base image as well as a context
   * to concoct the overlap image.
   * @param  {number}   precision   Scale for single hexagon-map image.
   * @param  {Array}    size        Width and height of base element.
   * @param  {function} pathGen     D3 path generator function.
   * @param  {Object}   geo         GeoJSON representing the object to tesselate.
   * @param  {number}   r           Hexagon radius.
   * @param  {number}   band        Extension of image (factor of r).
   * @return {Object}               The hex & geo image plus the context to use.
   */
  function getEdgeTools (precision, size, pathGen, geo, r, band) {
    // 1) Draw a hex with the correct radius at 0, 0.

    // Set up canvas and context.
    const w = Math.sqrt(3) * r * precision;
    const h = r * 2 * precision;
    const canvasHex = document.createElement('canvas');
    // const canvasHex = d3.select('body').append('canvas').node()
    canvasHex.width = w;
    canvasHex.height = h;
    const contextHex = canvasHex.getContext('2d');

    // Get the hexagon's corner points.
    const hexCorners = Array(7);
    for (let i = 0; i < 7; i++) {
      const corner = pointyHexCorner({
        x: 0,
        y: 0
      }, r * precision, i);
      hexCorners[i] = corner;
    }

    // Draw the hexagon.
    contextHex.translate(w / 2, h / 2);
    hexDraw(contextHex, hexCorners, 'red', 'fill');

    // 2) Draw the image.

    // Set up the image canvas and context.
    const [width, height] = size;
    const canvasImage = document.createElement('canvas');
    // const canvasImage = d3.select('body').append('canvas').node();
    canvasImage.width = width * precision;
    canvasImage.height = height * precision;
    const contextImage = canvasImage.getContext('2d');

    // Set the context for the path generator for use with Canvas.
    pathGen.context(contextImage);

    // Draw the image.
    const gridExtentStroke = band * r;
    contextImage.scale(precision, precision);
    contextImage.beginPath();
    pathGen(geo);
    contextImage.lineWidth = gridExtentStroke;
    contextImage.fillStyle = 'blue';
    contextImage.strokeStyle = 'blue';
    contextImage.stroke();
    contextImage.fill();

    // Reset the pathGenerators context.
    pathGen.context(null);

    // 3) Create context to combine images;

    const canvasMix = document.createElement('canvas');
    // const canvasMix = d3.select('body').append('canvas').node()
    canvasMix.width = w;
    canvasMix.height = h;
    const contextMix = canvasMix.getContext('2d');
    return {
      canvasHex,
      canvasImage,
      contextMix
    };
  }

  // Debug
  // import { pointyHexCorner, hexDraw } from './utils';

  /**
   * Calculates the cover for a single hexagon by
   * overlaying the map at the given position.
   * @param  {Object} edge      The datum representing the edge center.
   * @param  {Object} tools     The image and drawing tools
   *                            to create the overlap image.
   * @param  {number} precision The scaling factor for the image
   *                            at the given hex radius.
   * @param  {number} r         The hex radius. Required only for debugging.
   * @return {Object}           The cover updated egde center datum.
   */
  function getCover (edge, tools, precision, r) {
    const {
      canvasHex,
      canvasImage,
      contextMix
    } = tools;
    const w = canvasHex.width;
    const h = canvasHex.height;

    // // Debug ↓ --------------------------------------------------------------

    // // const r = 7;
    // const hexCorners = Array(7);
    // for (let i = 0; i < 7; i++) {
    //   const corner = pointyHexCorner({ x: 0, y: 0 }, r * precision, i);
    //   hexCorners[i] = corner;
    // }

    // const contextImage = canvasImage.getContext('2d');

    // // Centers.
    // contextImage.beginPath();
    //   contextImage.arc(edge.x, edge.y, 2, 0, 2*Math.PI)
    // contextImage.fillStyle = '#000'
    // contextImage.fill();

    // // Hexagons
    // hexDraw(contextImage, hexCorners, 'red', 'fill')

    // // Debug ↑ --------------------------------------------------------------

    // 1) Concoct the specific edge hexagon image and get the pixel data.

    // Draw hex image.
    contextMix.drawImage(canvasHex, 0, 0);

    // Set the composite type in preperation for the image overlap.
    contextMix.globalCompositeOperation = 'source-atop';

    // Draw Map at correct position.
    contextMix.drawImage(canvasImage, -edge.x * precision + w / 2, -edge.y * precision + h / 2);

    // Get the image data.
    const imageData = contextMix.getImageData(0, 0, w, h).data;

    // // Clear the canvas and reset the composite type in preperation
    // // for the next overlap (http://bit.do/ekDx4).
    // contextMix.clearRect(0,0,w,h);
    // contextMix.globalCompositeOperation = 'source-over';

    // 2) Calculate the image cover per edge hexagon.

    // Init area count variables.
    let hexArea = 0;
    let imgArea = 0;

    // Find filled pixel with some alpha (>=100)
    // and identify image part.
    for (let pixelIndex = 3; pixelIndex < imageData.length; pixelIndex += 4) {
      const alpha = imageData[pixelIndex];
      if (alpha < 100) {
        continue;
      } else {
        const red = imageData[pixelIndex - 3];
        const blue = imageData[pixelIndex - 1];
        red > blue ? hexArea++ : imgArea++;
      }
    }

    // Calculate cover and add to edge hexagon.
    const imgRatio = imgArea / (hexArea + imgArea);
    const updatedEdge = Object.assign({}, edge);
    updatedEdge.cover = imgRatio;

    // Clear the canvas and reset the composite type in preperation
    // for the next overlap (http://bit.do/ekDx4).
    contextMix.clearRect(0, 0, w, h);
    contextMix.globalCompositeOperation = 'source-over';
    return updatedEdge;
  }

  /**
   * Adds the updated cover value to each center datum.
   * @param  {Array}  centers   All center objects including the edge centers.
   * @param  {Array}  edges     Only the edge center objects.
   * @return {Array}            The updated center objects.
   */
  function addCover (centers, edges) {
    const centersUpdated = centers.slice(0);
    for (let i = 0; i < edges.length; i++) {
      const edge = edges[i];
      // Assuming the centers array id's are
      // consistent with the edge id's.
      centersUpdated[edge.id].cover = edge.cover;
    }
    return centersUpdated;
  }

  /**
   * Defines the data's latitude and longitude keys.
   * @param  {Array} lonLat   User defined array of geo keys.
   * @param  {Array} data     User defined data.
   * @return {Array}          Array of geo keys.
   */
  function checkGeoKeyNames(lonLat, data) {
    if (lonLat && lonLat.length === 2) return lonLat;
    const lonKey = Object.keys(data[0]).filter(key => {
      const low = key.toLowerCase();
      return low === 'longitude' || low === 'lon' || low === 'lng' || low === 'long' || low === 'lambda';
    });
    const latKey = Object.keys(data[0]).filter(key => {
      const low = key.toLowerCase();
      return low === 'latitude' || low === 'lat' || low === 'phi';
    });
    return [lonKey[0], latKey[0]];
  }

  /**
   * Process the user data to be structured for further use.
   * @param  {Array}    data          Array of user data objects.
   * @param  {function} projection    Geo projection.
   * @param  {Array}    variables     Optional. Array of variables the user
   *                                  would like to add to the layout.
   * @return {Array}                  Array of user's data points.
   */
  function prepUserData (data, projection, lonLat, variables) {
    // Return an empty array if the user hasn't passed down data.
    if (!data.length) return [];
    const geoKeys = checkGeoKeyNames(lonLat, data);
    return data.map(el => {
      const coords = projection([+el[geoKeys[0]], +el[geoKeys[1]]]);
      const obj = {};
      [obj.x, obj.y] = coords;
      if (variables && variables.length) {
        variables.forEach(varName => {
          obj[varName] = el[varName];
        });
      }
      return obj;
    });
  }

  /* eslint-disable no-param-reassign */

  /**
   * Bring each hexpoint into shape, by rolling up number of datapoints
   * per hexagon, add cover and setting apart original centers from
   * centers added by user-data.
   * @param  {Array} hexPoints        Array of arrays of grid and
   *                                  datapoints per hexagon.
   * @return {Array}                  Array of arrays of datapoints
   *                                  per hexagon plus additional props.
   */
  function rollupPoints (hexPoints) {
    for (let i = 0; i < hexPoints.length; i++) {
      // Cache current element and prep cover variable.
      const hexPoint = hexPoints[i];
      let cover;
      let gridpoint;

      // Remove grid points and cache cover.
      for (let j = 0; j < hexPoint.length; j++) {
        if (hexPoint[j].gridpoint === 1) {
          cover = hexPoint[j].cover;
          gridpoint = 1;
          hexPoint.splice(j, 1);
        }
      }

      // Augment with new properties.
      hexPoints[i].datapoints = hexPoints[i].length;
      hexPoints[i].cover = cover;
      hexPoints[i].gridpoint = gridpoint || 0;
    }
    return hexPoints;
  }

  /**
   * Calculates the cover weighted measures. Also assigns a
   * minimum cover proxy to each layout point without a cover.
   * Requried as some user data points can lie just outside the image.
   * @param  {Array}  points  Layout objects.
   * @param  {number} r       The hexagon's radius.
   * @return {Array}          Cover augmented layout objects.
   */
  function rollupDensity (points, r) {
    // Establish a minimum cover proxy: get a sorted array of cover values
    // for the quantile function. Only consider edges with cover < 1.
    const ascendingCover = points.filter(p => p.cover > 0 && p.cover < 1).map(d => d.cover).sort((a, b) => a - b);
    // Get the 10th percentile as the proxy.
    const quartileCover = quantile(ascendingCover, 0.1);

    // Get the hexagon's area in square pixel.
    const hexArea = 3 / 2 * Math.sqrt(3) * r ** 2;

    // Initialise extents.
    let maxPoints = 0;
    let maxPointsWt = 0;
    let maxDensity = 0;

    // Initialise the min values with the largest possible min value.
    let minPoints = points.length;
    let minPointsWt = points.length;
    let minDensity = points.length / hexArea;
    for (let i = 0; i < points.length; i++) {
      const point = points[i];

      // All layout points w/o cover will get assigned the cover proxy.
      // Note, only non-gridpoont datapoints will have no cover.
      if (!point.cover) {
        point.cover = quartileCover;
      }

      // Calculate the cover weighted measures.
      point.datapointsWt = point.datapoints * (1 / point.cover);
      point.pointDensity = point.datapoints / (hexArea * point.cover);

      // Update extents.
      maxPoints = Math.max(maxPoints, point.datapoints);
      maxPointsWt = Math.max(maxPointsWt, point.datapointsWt);
      maxDensity = Math.max(maxDensity, point.pointDensity);
      if (point.datapoints > 0) {
        minPoints = Math.min(minPoints, point.datapoints);
        minPointsWt = Math.min(minPointsWt, point.datapointsWt);
      }
      if (point.pointDensity > 0) {
        minDensity = Math.min(minDensity, point.pointDensity);
      }
    }
    const extentPoints = [minPoints, maxPoints];
    const extentPointsWeighted = [minPointsWt, maxPointsWt];
    const extentPointDensity = [minDensity, maxDensity];
    return {
      layout: points,
      extentPoints,
      extentPointsWeighted,
      extentPointDensity
    };
  }

  /**
   * Main hexgrid component.
   */
  function hexgrid () {
    // Init exposed.
    let extent;
    let geography;
    let projection;
    let pathGenerator;
    let hexRadius = 4;
    let hexRadiusUnit = null;
    let hexRadiusInUnits = null;
    let edgePrecision = 1;
    let gridExtend = 0;
    let geoKeys;

    /**
     * hexgrid function producing the layout.
     * @param  {Array} userData       Datapoints to visualise.
     *                                One datum represents one location.
     * @param  {Array} userVariables  Optional array of object keys to be
     *                                included in the final layout hex data.
     * @return {function/Object}      Augmented hexbin generator.
     */
    const hexgrid = function (userData, userVariables) {
      // Convert to pixel radius if provided in units.
      if (hexRadiusInUnits) {
        const conversion = convertUnitRadius(hexRadiusInUnits, hexRadiusUnit, projection, extent[1].map(d => d / 2));
        hexRadius = conversion.radiusPixel;
      }

      // Set hex radius to nearest full- or half-pixel.
      hexRadius = Math.round(hexRadius * 2) / 2;

      // Identify hexagons to draw.
      const hexbin = setHexGenerator(extent, hexRadius);
      const size = hexbin.size();
      const centers = hexbin.centers();
      const imageData = getImageData(size, pathGenerator, geography, hexRadius, 'fill', gridExtend);
      let imageCenters = getImageCenters(centers, imageData, size);

      // Identify edge hexagons and calculate image overlap ratio.
      const imageDataEdges = getImageData(size, pathGenerator, geography, hexRadius, 'stroke', gridExtend);
      const imageEdges = getEdgeCenters(imageCenters, imageDataEdges, size);
      const edgeTools = getEdgeTools(edgePrecision, size, pathGenerator, geography, hexRadius, gridExtend);
      const imageEdgesCover = imageEdges.map(d => getCover(d, edgeTools, edgePrecision));
      imageCenters = addCover(imageCenters, imageEdgesCover);

      // Prepare user data to augment layout.
      const userDataPrepped = prepUserData(userData, projection, geoKeys, userVariables);
      const mergedData = imageCenters.concat(userDataPrepped);
      const hexPoints = hexbin(mergedData);
      let hexData = rollupPoints(hexPoints);
      hexData = rollupDensity(hexData, hexRadius);

      // Augment hexbin generator.
      hexbin.grid = {};
      hexbin.grid.layout = hexData.layout;
      hexbin.grid.imageCenters = imageCenters;
      hexbin.grid.extentPoints = hexData.extentPoints;
      hexbin.grid.extentPointsWeighted = hexData.extentPointsWeighted;
      hexbin.grid.extentPointDensity = hexData.extentPointDensity;
      return hexbin;
    };

    // Exposed.
    hexgrid.extent = function (_) {
      return arguments.length ? (extent = expandExtent(_), hexgrid) : extent;
    };
    hexgrid.geography = function (_) {
      return arguments.length ? (geography = _, hexgrid) : geography;
    };
    hexgrid.projection = function (_) {
      return arguments.length ? (projection = _, hexgrid) : projection;
    };
    hexgrid.pathGenerator = function (_) {
      return arguments.length ? (pathGenerator = _, hexgrid) : pathGenerator;
    };
    hexgrid.hexRadius = function () {
      for (var _len = arguments.length, args = new Array(_len), _key = 0; _key < _len; _key++) {
        args[_key] = arguments[_key];
      }
      if (!args.length) {
        return hexRadiusUnit ? {
          radius: hexRadius,
          unit: hexRadiusUnit
        } : hexRadius;
      }
      if (args.length === 1) {
        return hexRadius = args[0], hexgrid;
      }
      if (args.length === 2) {
        [hexRadiusInUnits, hexRadiusUnit] = args;
        return hexgrid;
      }
      throw new Error('Please pass a numeric radius and optionally a string distance unit ("miles" or "kilometres") to `.hexradius()`');
    };
    hexgrid.edgePrecision = function (_) {
      return arguments.length ? (edgePrecision = convertToMin(_, 'Edge precision', 0.3), hexgrid) : edgePrecision;
    };
    hexgrid.gridExtend = function (_) {
      return arguments.length ? (gridExtend = convertToMin(_, 'Edge band', 0), hexgrid) : gridExtend;
    };
    hexgrid.geoKeys = function (_) {
      return arguments.length ? (geoKeys = _, hexgrid) : geoKeys;
    };
    return hexgrid;
  }

  /**
   * Produce an array or arrays with all polygonal boundary points.
   * @param  {Object}   geo           The GeoJSON FeatureCollection.
   * @param  {function} projection    The D3 projection function.
   * @return {Array}                  Array of arrays holding the boundary points
   *                                  for each area.
   */
  function getBoundaryPoints (geo, projection) {
    let boundaryPoints = [];
    let collection;

    // 1) Try for geometry type and get their contents.

    try {
      if (geo.type === 'FeatureCollection') {
        collection = geo.features;
      } else if (geo.type === 'GeometryCollection') {
        collection = geo.geometries;
      } else {
        throw new Error('Geometry type not supported. Please feed me a "FeatureCollection" or a "GeometryCollection".');
      }
    } catch (err) {
      throw new Error(err);
    }

    // 2) Retrieve the boundary points.

    for (let i = 0; i < collection.length; i++) {
      // Crack open the geometry to get the coordinate holder object.
      const geom = geo.type === 'FeatureCollection' ? geo.features[i].geometry : geo.geometries[i];

      // Different ways to access coordinates in a FeatureCollection:

      // Polygons: coordinates[Array[coordinates]]
      if (geom && geom.type === 'Polygon') {
        // Correcting for longitudes +180°.
        const polygon = geom.coordinates[0].map(coord => projection(coord[0] > 180 ? [180, coord[1]] : coord));
        boundaryPoints.push(polygon);

        // MultiPolygons: coordinates[Polygons[Array[[coordinates]]]]
      } else if (geom && geom.type === 'MultiPolygon') {
        // Correcting for longitudes +180°.
        const polygons = geom.coordinates.map(multi => multi[0].map(coord => projection(coord[0] > 180 ? [180, coord[1]] : coord)));
        boundaryPoints = boundaryPoints.concat(polygons);
      } else {
        continue;
      }
    }
    return boundaryPoints;
  }

  function polygonContains(polygon, point) {
    var n = polygon.length,
        p = polygon[n - 1],
        x = point[0], y = point[1],
        x0 = p[0], y0 = p[1],
        x1, y1,
        inside = false;

    for (var i = 0; i < n; ++i) {
      p = polygon[i], x1 = p[0], y1 = p[1];
      if (((y1 > y) !== (y0 > y)) && (x < (x0 - x1) * (y - y1) / (y0 - y1) + x1)) inside = !inside;
      x0 = x1, y0 = y1;
    }

    return inside;
  }

  /**
   * Produce an array or arrays with all points within a polygonial area/feature.
   * @param  {Array} gridPoints         All grid points.
   * @param  {Array} boundaryPoints     Array of arrays, one for each area,
   *                                    holding the area's boundary points.
   * @return {Array}                    Array of grid points within each area.
   *                                    Sorted ascendingly by x and y.
   */
  function getPolygonPoints (gridPoints, boundaryPoints) {
    return boundaryPoints.reduce((result, boundary) => {
      const areaPoints = gridPoints.filter(point => polygonContains(boundary, [point.x, point.y]));
      return result.concat(areaPoints);
    }, []).sort((a, b) => a.x - b.x || a.y - b.y);
  }

  exports.geoPolygon = getBoundaryPoints;
  exports.hexgrid = hexgrid;
  exports.polygonPoints = getPolygonPoints;

  Object.defineProperty(exports, '__esModule', { value: true });

}));
